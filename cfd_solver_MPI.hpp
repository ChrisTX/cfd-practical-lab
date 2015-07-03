#ifndef CFD_SOLVER_HPP
#define CFD_SOLVER_HPP

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include <functional>
#include <limits>
#include <utility>
#include <map>
#include <set>
#include <vector>

#include "cfd_vtk_library.hpp"
#include "fdm_library.hpp"
#include "linalg_vector2d.hpp"

#include <iostream>

#include <mpi.h>
#include "mpi_helper.hpp"

namespace CFD {

	template<
		typename T,
		template<typename...> class VectorT = LinearAlgebra::OffsetVector2D
	> class CFDSolverMPI {
	public:
		using VectorType = VectorT<T>;
		using size_type = typename VectorType::size_type;

		static_assert(std::is_floating_point<T>(), "T has to be a floating point type");
		static_assert(std::numeric_limits<T>::has_signaling_NaN, "We require a signaling NaN");

		using input_func_type = std::function<std::pair<T, T>(T x, T y)>;

		enum class BoundaryCondition_t {
			noslip,
			freeslip,
			outflow,
			inflow
		};

		enum class BoundaryEdge_t {
			left,
			right,
			top,
			bottom
		};

		struct numeric_parameters {
			T γ;
			T ω;
			T ε;
			size_type itermax;
		};

		struct problem_data {
			T Re;
			T gx;
			T gy;
		};

		struct time_parameters {
			T τ; ///< safety factor for time step size control
			T t; ///< current time value
			T δt; ///< if τ < 0, this will be used as fixed step size
			T t_end; ///< final time
			T t_vis; ///< time step for visualization
		};

		struct problem_geometry {
			T xlength;
			T ylength;
			size_type imax;
			size_type jmax;
		};

		struct MPI_information {
			int rank;
			MPI_Comm commhandle;
			size_type idimmax;
			size_type jdimmax;
		};

		CFDSolverMPI(MPI_information mpiinfo, problem_data pbdata, problem_geometry pbgeom, time_parameters timeparam, numeric_parameters numparam, const std::function<T(T,T)>& u0_init, const std::function<T(T,T)>& v0_init, const std::function<T(T,T)>& p0_init)
		: m_MPIInfo(std::move(mpiinfo)), m_ProblemData(std::move(pbdata)), m_ProblemGeometry(std::move(pbgeom)), m_TimeParameters(std::move(timeparam)), m_NumericParameters(std::move(numparam))
		{
			DetermineSpecificMPIData();
			const size_type sizeI = m_MPISpecData.ir - m_MPISpecData.il + 1;
			const size_type sizeJ = m_MPISpecData.jt - m_MPISpecData.jb + 1;
			m_u = VectorType(sizeI + 2, sizeJ + 2, m_MPISpecData.il - 1, m_MPISpecData.jb - 1);
			m_v = VectorType(sizeI + 2, sizeJ + 2, m_MPISpecData.il - 1, m_MPISpecData.jb - 1);
			m_p = VectorType(sizeI + 2, sizeJ + 2, m_MPISpecData.il - 1, m_MPISpecData.jb - 1);

			for(size_type i = m_MPISpecData.il - 1; i <= m_MPISpecData.ir + 1; ++i) {
				for(size_type j = m_MPISpecData.jb - 1; j <= m_MPISpecData.jt + 1; ++j) {
					if(i && i <= m_ProblemGeometry.imax && j && j <= m_ProblemGeometry.jmax) {
						const T ux = δx() * (i - 1);
						const T uy = ( δy() * (j - 1) ) / 2;
						m_u(i, j) = u0_init(ux, uy);

						const T vx = ( δx() * (i - 1) ) / 2;
						const T vy = δy() * (j - 1);
						m_v(i, j) = v0_init(vx, vy);

						m_p(i, j) = p0_init(vx, uy);
					}
				}
			}
		}

		void Solve()
		{
			T t_since_last_print = 0;
			VTK::TimedVTKRectilinearGridPrinter<T> CFD_VTKPrinters("output-", 0u, "CFD solution", VTK::VTKFileFormat_t::ASCII, static_cast<VTK::VTKIndex_t>(m_ProblemGeometry.imax), static_cast<VTK::VTKIndex_t>(m_ProblemGeometry.jmax), δx(), δy() );
			while(m_TimeParameters.t < m_TimeParameters.t_end) {
				ApplyBoundaryConditions();
				SolveLoop();
				t_since_last_print += δt;
				if(t_since_last_print >= m_TimeParameters.t_vis) {
					t_since_last_print = 0;

					if(m_MPIInfo.rank == 0) {

						auto u = VectorType(m_ProblemGeometry.imax + 2, m_ProblemGeometry.jmax + 2, 0u, 0u, std::numeric_limits<T>::signaling_NaN());
						auto v = VectorType(m_ProblemGeometry.imax + 2, m_ProblemGeometry.jmax + 2, 0u, 0u, std::numeric_limits<T>::signaling_NaN());
						auto p = VectorType(m_ProblemGeometry.imax + 2, m_ProblemGeometry.jmax + 2, 0u, 0u, std::numeric_limits<T>::signaling_NaN());

						for(int i = 0; i < static_cast<int>( m_MPIInfo.idimmax ); ++i) {
							for(int j = 0; j < static_cast<int>( m_MPIInfo.jdimmax ); ++j) {
								const int target_coords[2u] = { i, j };
								int target_rank;
								if( MPI_Cart_rank(m_MPIInfo.commhandle, target_coords, &target_rank) != MPI_SUCCESS )
									throw std::runtime_error{"MPI_Cart_rank failed!"};

								size_type sizeI = m_ProblemGeometry.imax/m_MPIInfo.idimmax + 2;
								size_type sizeJ = m_ProblemGeometry.jmax/m_MPIInfo.jdimmax + 2;
								if( i + 1 == m_MPIInfo.idimmax )
									sizeI += m_ProblemGeometry.imax % m_MPIInfo.idimmax;
								if( j + 1 == m_MPIInfo.jdimmax )
									sizeJ += m_ProblemGeometry.jmax % m_MPIInfo.jdimmax;

								size_type ioff = m_ProblemGeometry.imax/m_MPIInfo.idimmax * i;
								size_type joff = m_ProblemGeometry.jmax/m_MPIInfo.jdimmax * j;

								if(target_rank == 0) {
									// We're assuming row-major storage for TVector's underlying data here. A bit dangerous as we can't ensure this.
									
									for(size_type k = 1u; k + 1 < sizeJ; ++k) {
										for(size_type l = 1u; l + 1 < sizeI; ++l) {
											assert(std::isnan(u(ioff + l, joff + k)));
											assert(std::isnan(v(ioff + l, joff + k)));
											assert(std::isnan(p(ioff + l, joff + k)));
											u(ioff + l, joff + k) = m_u.data()[l + k * sizeI];
											v(ioff + l, joff + k) = m_v.data()[l + k * sizeI];
											p(ioff + l, joff + k) = m_p.data()[l + k * sizeI];
										}
									}
								}
								else {
									size_type totalsize = sizeI * sizeJ;
									std::vector<T> local_u(totalsize);
									std::vector<T> local_v(totalsize);
									std::vector<T> local_p(totalsize);
									MPI_Request mpi_reqs[3];
									if( MPI_Irecv(local_u.data(), static_cast<int>( totalsize ), MPIHelper::MPIDatatype_t<T>::value(), target_rank, static_cast<int>( MPITags_t::VTKU ), m_MPIInfo.commhandle, mpi_reqs) != MPI_SUCCESS )
										throw std::runtime_error{"MPI_Irecv failed!"};
									if( MPI_Irecv(local_v.data(), static_cast<int>( totalsize ), MPIHelper::MPIDatatype_t<T>::value(), target_rank, static_cast<int>( MPITags_t::VTKV ), m_MPIInfo.commhandle, mpi_reqs + 1) != MPI_SUCCESS )
										throw std::runtime_error{"MPI_Irecv failed!"};
									if( MPI_Irecv(local_p.data(), static_cast<int>( totalsize ), MPIHelper::MPIDatatype_t<T>::value(), target_rank, static_cast<int>( MPITags_t::VTKP ), m_MPIInfo.commhandle, mpi_reqs + 2) != MPI_SUCCESS )
										throw std::runtime_error{"MPI_Irecv failed!"};

									if( MPI_Waitall(3, mpi_reqs, MPI_STATUSES_IGNORE) != MPI_SUCCESS )
										throw std::runtime_error{"MPI_Waitall failed!"};

									// We're assuming row-major storage for TVector's underlying data here. A bit dangerous as we can't ensure this.
									for(size_type k = 1u; k + 1 < sizeJ; ++k) {
										for(size_type l = 1u; l + 1 < sizeI; ++l) {
											assert(std::isnan(u(ioff + l, joff + k)));
											assert(std::isnan(v(ioff + l, joff + k)));
											assert(std::isnan(p(ioff + l, joff + k)));
											u(ioff + l, joff + k) = local_u[l + k * sizeI];
											v(ioff + l, joff + k) = local_v[l + k * sizeI];
											p(ioff + l, joff + k) = local_p[l + k * sizeI];
										}
									}
								}
							}
						}

						auto& vtkprt = CFD_VTKPrinters.GetNextPrinter();
						vtkprt.UsePointData();
						vtkprt.AddScalars("Pressure", p);
						vtkprt.AddVectors("Velocity", u, v);
					}
					else {
						MPI_Request mpi_reqs[3];
						if( MPI_Isend(m_u.data(), static_cast<int>( m_u.size() ), MPIHelper::MPIDatatype_t<T>::value(), 0, static_cast<int>( MPITags_t::VTKU ), m_MPIInfo.commhandle, mpi_reqs) != MPI_SUCCESS )
							throw std::runtime_error{"MPI_Isend failed!"};
						if( MPI_Isend(m_v.data(), static_cast<int>( m_v.size() ), MPIHelper::MPIDatatype_t<T>::value(), 0, static_cast<int>( MPITags_t::VTKV ), m_MPIInfo.commhandle, mpi_reqs + 1) != MPI_SUCCESS )
							throw std::runtime_error{"MPI_Isend failed!"};
						if( MPI_Isend(m_p.data(), static_cast<int>( m_p.size() ), MPIHelper::MPIDatatype_t<T>::value(), 0, static_cast<int>( MPITags_t::VTKP ), m_MPIInfo.commhandle, mpi_reqs + 2) != MPI_SUCCESS )
							throw std::runtime_error{"MPI_Isend failed!"};

						if( MPI_Waitall(3, mpi_reqs, MPI_STATUSES_IGNORE) != MPI_SUCCESS )
							throw std::runtime_error{"MPI_Waitall failed!"};
					}
				}
			}
		}

		void SetBoundaryCondition(BoundaryEdge_t bc_edge, BoundaryCondition_t bc_type)
		{
			m_BCs[bc_edge] = std::make_pair(bc_type, nullptr);
		}

		void SetBoundaryCondition(BoundaryEdge_t bc_edge, BoundaryCondition_t bc_type, const input_func_type& bc_func)
		{
			if(bc_type != BoundaryCondition_t::inflow)
				throw std::invalid_argument("Inflow function on other boundary type given!");

			m_BCs[bc_edge] = std::make_pair(bc_type, bc_func);
		}

	private:
		MPI_information m_MPIInfo;

		struct MPI_specific_data
		{
			size_type ir;
			size_type il;
			size_type jt;
			size_type jb;
			size_type ipos;
			size_type jpos;
		};

		MPI_specific_data m_MPISpecData;

		T δt; ///< time step size
		inline T δx() const { return m_ProblemGeometry.xlength/m_ProblemGeometry.imax; } 
		inline T δy() const { return m_ProblemGeometry.ylength/m_ProblemGeometry.jmax; }


		void DetermineSpecificMPIData()
		{
			int my_coords[2];
			if(MPI_Cart_coords(m_MPIInfo.commhandle, m_MPIInfo.rank, 2, my_coords) != MPI_SUCCESS)
				throw std::runtime_error{"MPI_Cart_coords failed"};

			m_MPISpecData.ipos = static_cast<size_type>( my_coords[0] );
			m_MPISpecData.jpos = static_cast<size_type>( my_coords[1] );

			assert(m_MPISpecData.ipos < m_MPIInfo.idimmax && m_MPISpecData.jpos < m_MPIInfo.jdimmax);

			// The following is intended to be integer division!
			m_MPISpecData.ir = m_ProblemGeometry.imax/m_MPIInfo.idimmax * (m_MPISpecData.ipos + 1);
			m_MPISpecData.il = m_ProblemGeometry.imax/m_MPIInfo.idimmax * m_MPISpecData.ipos + 1;
			m_MPISpecData.jt = m_ProblemGeometry.jmax/m_MPIInfo.jdimmax * (m_MPISpecData.jpos + 1);
			m_MPISpecData.jb = m_ProblemGeometry.jmax/m_MPIInfo.jdimmax * m_MPISpecData.jpos + 1;

			if(m_ProblemGeometry.imax % m_MPIInfo.idimmax && m_MPISpecData.ipos + 1 == m_MPIInfo.idimmax) {
				assert(m_ProblemGeometry.imax == m_MPISpecData.ir + m_ProblemGeometry.imax % m_MPIInfo.idimmax);
				m_MPISpecData.ir = m_ProblemGeometry.imax;
			}

			if(m_ProblemGeometry.jmax % m_MPIInfo.jdimmax && m_MPISpecData.jpos + 1 == m_MPIInfo.jdimmax) {
				assert(m_ProblemGeometry.jmax == m_MPISpecData.jt + m_ProblemGeometry.jmax % m_MPIInfo.jdimmax);
				m_MPISpecData.jt = m_ProblemGeometry.jmax;
			}

			assert(m_MPISpecData.il && m_MPISpecData.jb && (m_MPISpecData.ir <= m_ProblemGeometry.imax) && (m_MPISpecData.jt <= m_ProblemGeometry.jmax));
		}

		bool IsRootProcess() const
		{
			return (m_MPIInfo.rank == 0);
		}

		void Updateδt()
		{
			if(m_TimeParameters.τ < 0) {
				δt = m_TimeParameters.δt;
				return;
			}

			T local_max[2] = {T{0}, T{0}};
			for(size_type i = m_MPISpecData.il; i <= m_MPISpecData.ir; ++i) {
				for(size_type j = m_MPISpecData.jb; j <= m_MPISpecData.jt; ++j) {
					local_max[0] = std::max( local_max[0], std::abs(m_u(i, j)) );
					local_max[1] = std::max( local_max[1], std::abs(m_v(i, j)) );
				}
			}

			T global_max[2];
			if(MPI_Allreduce(&local_max, &global_max, 2, MPIHelper::MPIDatatype_t<T>::value(), MPI_MAX, m_MPIInfo.commhandle) != MPI_SUCCESS)
				throw std::runtime_error{"MPI_Allreduce failed!"};

			δt = m_TimeParameters.τ * std::min( {
				( m_ProblemData.Re / T{2} ) * std::pow( T{1} / std::pow(δx(), 2) + T{1} / std::pow(δy(), 2), -1),
				δx() / global_max[0],
				δy() / global_max[1]
			} );
		}

		void CalculateFG()
		{
			m_F.resize(m_MPISpecData.ir - m_MPISpecData.il + 2, m_MPISpecData.jt - m_MPISpecData.jb + 2, std::numeric_limits<T>::signaling_NaN());
			m_F.SetOffset(m_MPISpecData.il - 1, m_MPISpecData.jb - 1);
			m_G.resize(m_MPISpecData.ir - m_MPISpecData.il + 2, m_MPISpecData.jt - m_MPISpecData.jb + 2, std::numeric_limits<T>::signaling_NaN());
			m_G.SetOffset(m_MPISpecData.il - 1, m_MPISpecData.jb - 1);

			const auto FDM = FDM::FDMVectors2D<VectorType>(m_NumericParameters.γ, δx(), δy(), m_u, m_v);

			for(size_type i = (m_MPISpecData.il == 1 ? 0u : m_MPISpecData.il) ; i <= m_MPISpecData.ir; ++i) {
				for(size_type j = (m_MPISpecData.jb == 1 ? 0u : m_MPISpecData.jb); j <= m_MPISpecData.jt; ++j) {
					if(j != m_MPISpecData.jb - 1) {
						m_F(i, j) = m_u(i, j);

						if(i && i < m_ProblemGeometry.imax)
							m_F(i, j) += δt * ( ( T{1} / m_ProblemData.Re ) * ( FDM.d2u_dx2(i, j) + FDM.d2u_dy2(i, j) ) - FDM.du2_dx(i, j) - FDM.duv_dy(i, j) + m_ProblemData.gx );
						
						assert(std::isfinite(m_F(i, j)));
					}

					if(i != m_MPISpecData.il - 1) {
						m_G(i, j) = m_v(i, j);

						if(j && j < m_ProblemGeometry.jmax)
							m_G(i, j) += δt * ( ( T{1} / m_ProblemData.Re ) * ( FDM.d2v_dx2(i, j) + FDM.d2v_dy2(i, j) ) - FDM.duv_dx(i, j) - FDM.dv2_dy(i, j) + m_ProblemData.gy );
						
						assert(std::isfinite(m_G(i, j)));
					}
				}
			}

			ExchangeFGData();
		}

		VectorType CalculateRHS() const
		{
			VectorType RHS(m_MPISpecData.ir - m_MPISpecData.il + 1, m_MPISpecData.jt - m_MPISpecData.jb + 1, m_MPISpecData.il, m_MPISpecData.jb, std::numeric_limits<T>::signaling_NaN());

			for(size_type i = m_MPISpecData.il; i <= m_MPISpecData.ir; ++i) {
				for(size_type j = m_MPISpecData.jb; j <= m_MPISpecData.jt; ++j) {
					RHS(i, j) = (T{1} / δt) * ( ( m_F(i, j) - m_F(i - 1, j) )/δx() + ( m_G(i, j) - m_G(i, j - 1) )/δy() );
					assert(std::isfinite(RHS(i, j)));
				}
			}

			return RHS;
		}

		T CalculateResidual(const VectorType& RHS) const
		{
			T local_residual{0};
			for(size_type i = m_MPISpecData.il; i <= m_MPISpecData.ir; ++i) {
				for(size_type j = m_MPISpecData.jb; j <= m_MPISpecData.jt; ++j) {
					local_residual += std::pow( ( m_p(i+1, j) - 2 * m_p(i, j) + m_p(i-1, j) )/( std::pow(δx(), 2) )
										+ ( m_p(i, j+1) - 2 * m_p(i, j) + m_p(i, j-1) )/( std::pow(δy(), 2) )
										- RHS(i, j), 2);
				}
			}
			assert(std::isfinite(local_residual));

			T residual;
			if(MPI_Allreduce(&local_residual, &residual, 1, MPIHelper::MPIDatatype_t<T>::value(), MPI_SUM, m_MPIInfo.commhandle) != MPI_SUCCESS)
				throw std::runtime_error{"MPI_Allreduce failed!"};

			assert(std::isfinite(residual));

			residual = std::sqrt( residual/(m_ProblemGeometry.imax * m_ProblemGeometry.jmax) );
			std::cout << residual << std::endl;
			return residual;
		}

		void SolveLoop()
		{
			Updateδt();
			CalculateFG();
			const auto RHS = CalculateRHS();

			for(size_type it_count = 0u; it_count < m_NumericParameters.itermax; ++it_count) {
				ApplyPressureConditions();

				for(size_type i = m_MPISpecData.il; i <= m_MPISpecData.ir; ++i) {
					for(size_type j = m_MPISpecData.jb; j <= m_MPISpecData.jt; ++j) {
						m_p(i,j) = (1 - m_NumericParameters.ω) * m_p(i, j) + ( m_NumericParameters.ω / ( T{2} / std::pow(δx(), 2) + T{2} / std::pow(δy(), 2) ) ) *
								 ( ( m_p(i+1,j) + m_p(i-1, j) ) / std::pow(δx(), 2) + ( m_p(i,j+1) + m_p(i, j-1) ) / std::pow(δy(), 2) - RHS(i, j) );
						assert(std::isfinite(m_p(i,j)));
					}
				}

				if(CalculateResidual(RHS) < m_NumericParameters.ε)
					break;
			}

			for(size_type i = m_MPISpecData.il; i <= std::min(m_MPISpecData.ir, m_ProblemGeometry.imax); ++i)
				for(size_type j = m_MPISpecData.jb; j <= m_MPISpecData.jt; ++j)
					m_u(i,j) = m_F(i,j) - δt/δx() * (m_p(i+1,j) - m_p(i,j));

			for(size_type i = m_MPISpecData.il; i <= m_MPISpecData.ir; ++i)
				for(size_type j = m_MPISpecData.jb; j <= std::min(m_MPISpecData.jt, m_ProblemGeometry.jmax); ++j)
					m_v(i,j) = m_G(i,j) - δt/δy() * (m_p(i,j+1) - m_p(i,j));

			ExchangeVelocityData();

			m_TimeParameters.t += δt;
		}

		std::map<BoundaryEdge_t, std::pair<BoundaryCondition_t, input_func_type>> m_BCs;

		void ApplyBoundaryConditions()
		{
			if(m_MPISpecData.jb == 1) {
				for(std::size_t i = m_MPISpecData.il; i <= m_MPISpecData.ir; ++i) {
					T u_value = std::numeric_limits<T>::signaling_NaN(), v_value = u_value;
					const auto bcpair = m_BCs.at(BoundaryEdge_t::bottom);
					switch(bcpair.first) {
						case BoundaryCondition_t::noslip:
							u_value = -m_u(i, 1);
							v_value = 0;
							break;

						case BoundaryCondition_t::freeslip:
							u_value = m_u(i, 1);
							v_value = 0;
							break;

						case BoundaryCondition_t::outflow:
							u_value = m_u(i, 1);
							v_value = m_v(i, 1);
							break;

						case BoundaryCondition_t::inflow:
							const auto input_fn_val = bcpair.second(i * δx(), 0);
							u_value = 2 * input_fn_val.first - m_u(i, 1);
							v_value = input_fn_val.second;
					}

					m_u(i, 0) = u_value;
					m_v(i, 0) = v_value;
				}
			}

			if(m_MPISpecData.jt == m_ProblemGeometry.jmax) {
				for(std::size_t i = m_MPISpecData.il; i <= m_MPISpecData.ir; ++i) {
					T u_value = std::numeric_limits<T>::signaling_NaN(), v_value = u_value;
					const auto bcpair = m_BCs.at(BoundaryEdge_t::top);
					switch(bcpair.first) {
						case BoundaryCondition_t::noslip:
							u_value = -m_u(i, m_ProblemGeometry.jmax);
							v_value = 0;
							break;

						case BoundaryCondition_t::freeslip:
							u_value = m_u(i, m_ProblemGeometry.jmax);
							v_value = 0;
							break;
							
						case BoundaryCondition_t::outflow:
							u_value = m_u(i, m_ProblemGeometry.jmax);
							v_value = m_v(i, m_ProblemGeometry.jmax - 1);
							break;
							
						case BoundaryCondition_t::inflow:
							const auto input_fn_val = bcpair.second(i * δx(), m_ProblemGeometry.ylength);
							u_value = 2 * input_fn_val.first - m_u(i, m_ProblemGeometry.jmax);
							v_value = input_fn_val.second;
					}

					m_u(i, m_ProblemGeometry.jmax + 1) = u_value;
					m_v(i, m_ProblemGeometry.jmax) = v_value;
				}
			}

			if(m_MPISpecData.il == 1) {
				for(std::size_t j = m_MPISpecData.jb; j <= m_MPISpecData.jt; ++j) {
					T u_value = std::numeric_limits<T>::signaling_NaN(), v_value = u_value;
					const auto bcpair = m_BCs.at(BoundaryEdge_t::left);
					switch(bcpair.first) {
						case BoundaryCondition_t::noslip:
							u_value = 0;
							v_value = -m_v(1, j);
							break;

						case BoundaryCondition_t::freeslip:
							u_value = 0;
							v_value = m_v(1, j);
							break;
							
						case BoundaryCondition_t::outflow:
							u_value = m_u(1, j);
							v_value = m_v(1, j);
							break;
							
						case BoundaryCondition_t::inflow:
							const auto input_fn_val = bcpair.second(0, j * δy());
							u_value = input_fn_val.first;
							v_value = 2 * input_fn_val.second - m_v(1, j);
					}

					m_u(0, j) = u_value;
					m_v(0, j) = v_value;
				}
			}

			if(m_MPISpecData.ir == m_ProblemGeometry.imax) {
				for(std::size_t j = m_MPISpecData.jb; j <= m_MPISpecData.jt; ++j) {
					T u_value = std::numeric_limits<T>::signaling_NaN(), v_value = u_value;
					const auto bcpair = m_BCs.at(BoundaryEdge_t::right);
					switch(bcpair.first) {
						case BoundaryCondition_t::noslip:
							u_value = 0;
							v_value = -m_v(m_ProblemGeometry.imax, j);
							break;

						case BoundaryCondition_t::freeslip:
							u_value = 0;
							v_value = m_v(m_ProblemGeometry.imax, j);
							break;
							
						case BoundaryCondition_t::outflow:
							u_value = m_u(m_ProblemGeometry.imax - 1, j);
							v_value = m_v(m_ProblemGeometry.imax, j);
							break;
							
						case BoundaryCondition_t::inflow:
							const auto input_fn_val = bcpair.second(m_ProblemGeometry.xlength, j * δy());
							u_value = input_fn_val.first;
							v_value = 2 * input_fn_val.second - m_v(m_ProblemGeometry.imax, j);
					}

					m_u(m_ProblemGeometry.imax, j) = u_value;
					m_v(m_ProblemGeometry.imax + 1, j) = v_value;
				}
			}
		}

		enum class MPITags_t {
			PressureRight,
			PressureLeft,
			PressureTop,
			PressureBottom,
			UVRight,
			UVLeft,
			UVTop,
			UVBottom,
			VTKU,
			VTKV,
			VTKP,
			F,
			G
		};

		void ApplyPressureConditions()
		{
			std::vector<MPI_Request> used_requests;
			std::vector<std::function<void()>> packaged_responses;
			std::vector<std::unique_ptr<T[]>> buffervector;

			if(m_MPISpecData.il == 1u) {
				for(size_type j = m_MPISpecData.jb; j <= m_MPISpecData.jt; ++j)
					m_p(0u, j) = m_p(1u, j);
			}
			else {
				int rank_source, rank_dest;
				if(MPI_Cart_shift(m_MPIInfo.commhandle, 0, -1, &rank_source, &rank_dest) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Cart_shift failed!"};

				const size_type transfer_size = m_MPISpecData.jt - m_MPISpecData.jb + 1;
				auto send_buffer = std::make_unique<T[]>(transfer_size);
				auto recv_buffer = std::make_unique<T[]>(transfer_size);
				for(size_type j = 0u; j < transfer_size; ++j)
					send_buffer.get()[j] = m_p(m_MPISpecData.il, j + m_MPISpecData.jb);
				
				MPI_Request send_request, recv_request;
				if(MPI_Isend(send_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::PressureRight ), m_MPIInfo.commhandle, &send_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Isend failed!"};
				if(MPI_Irecv(recv_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::PressureLeft ), m_MPIInfo.commhandle, &recv_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Irecv failed!"};

				used_requests.emplace_back(std::move(send_request));
				used_requests.emplace_back(std::move(recv_request));
				packaged_responses.emplace_back([recv_buf = recv_buffer.get(), this, transfer_size]() { 
					for(size_type j = 0u; j < transfer_size; ++j)
						this->m_p(m_MPISpecData.il - 1, j + m_MPISpecData.jb) = recv_buf[j];
				});
				buffervector.emplace_back(std::move(send_buffer));
				buffervector.emplace_back(std::move(recv_buffer));
			}

			if(m_MPISpecData.ir == m_ProblemGeometry.imax) {
				for(size_type j = m_MPISpecData.jb; j <= m_MPISpecData.jt; ++j)
					m_p(m_ProblemGeometry.imax + 1, j) = m_p(m_ProblemGeometry.imax, j);
			}
			else {
				int rank_source, rank_dest;
				if(MPI_Cart_shift(m_MPIInfo.commhandle, 0, 1, &rank_source, &rank_dest) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Cart_shift failed!"};

				const size_type transfer_size = m_MPISpecData.jt - m_MPISpecData.jb + 1;
				auto send_buffer = std::make_unique<T[]>(transfer_size);
				auto recv_buffer = std::make_unique<T[]>(transfer_size);
				for(size_type j = 0u; j < transfer_size; ++j)
					send_buffer.get()[j] = m_p(m_MPISpecData.ir, j + m_MPISpecData.jb);
				
				MPI_Request send_request, recv_request;
				if(MPI_Isend(send_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::PressureLeft ), m_MPIInfo.commhandle, &send_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Isend failed!"};
				if(MPI_Irecv(recv_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::PressureRight ), m_MPIInfo.commhandle, &recv_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Irecv failed!"};

				used_requests.emplace_back(std::move(send_request));
				used_requests.emplace_back(std::move(recv_request));
				packaged_responses.emplace_back([recv_buf = recv_buffer.get(), this, transfer_size]() { 
					for(size_type j = 0u; j < transfer_size; ++j)
						this->m_p(m_MPISpecData.ir + 1, j + m_MPISpecData.jb) = recv_buf[j];
				});
				buffervector.emplace_back(std::move(send_buffer));
				buffervector.emplace_back(std::move(recv_buffer));
			}

			if(m_MPISpecData.jb == 1u) {
				for(size_type i = m_MPISpecData.il; i <= m_MPISpecData.ir; ++i)
					m_p(i, 0u) = m_p(i, 1u);
			}
			else {
				int rank_source, rank_dest;
				if(MPI_Cart_shift(m_MPIInfo.commhandle, 1, -1, &rank_source, &rank_dest) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Cart_shift failed!"};

				const size_type transfer_size = m_MPISpecData.ir - m_MPISpecData.il + 1;
				auto send_buffer = std::make_unique<T[]>(transfer_size);
				auto recv_buffer = std::make_unique<T[]>(transfer_size);
				for(size_type i = 0u; i < transfer_size; ++i)
					send_buffer.get()[i] = m_p(i + m_MPISpecData.il, m_MPISpecData.jb);
				
				MPI_Request send_request, recv_request;
				if(MPI_Isend(send_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::PressureTop ), m_MPIInfo.commhandle, &send_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Isend failed!"};
				if(MPI_Irecv(recv_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::PressureBottom ), m_MPIInfo.commhandle, &recv_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Irecv failed!"};

				used_requests.emplace_back(std::move(send_request));
				used_requests.emplace_back(std::move(recv_request));
				packaged_responses.emplace_back([recv_buf = recv_buffer.get(), this, transfer_size]() { 
					for(size_type i = 0u; i < transfer_size; ++i)
						this->m_p(i + m_MPISpecData.il, m_MPISpecData.jb - 1) = recv_buf[i];
				});
				buffervector.emplace_back(std::move(send_buffer));
				buffervector.emplace_back(std::move(recv_buffer));
			}

			if(m_MPISpecData.jt == m_ProblemGeometry.jmax) {
				for(size_type i = m_MPISpecData.il; i <= m_MPISpecData.ir; ++i)
					m_p(i, m_ProblemGeometry.jmax + 1) = m_p(i, m_ProblemGeometry.jmax);
			}
			else {
				int rank_source, rank_dest;
				if(MPI_Cart_shift(m_MPIInfo.commhandle, 1, 1, &rank_source, &rank_dest) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Cart_shift failed!"};

				const size_type transfer_size = m_MPISpecData.ir - m_MPISpecData.il + 1;
				auto send_buffer = std::make_unique<T[]>(transfer_size);
				auto recv_buffer = std::make_unique<T[]>(transfer_size);
				for(size_type i = 0u; i < transfer_size; ++i)
					send_buffer.get()[i] = m_p(i + m_MPISpecData.il, m_MPISpecData.jt);
				
				MPI_Request send_request, recv_request;
				if(MPI_Isend(send_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::PressureBottom ), m_MPIInfo.commhandle, &send_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Isend failed!"};
				if(MPI_Irecv(recv_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::PressureTop ), m_MPIInfo.commhandle, &recv_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Irecv failed!"};

				used_requests.emplace_back(std::move(send_request));
				used_requests.emplace_back(std::move(recv_request));
				packaged_responses.emplace_back([recv_buf = recv_buffer.get(), this, transfer_size]() { 
					for(size_type i = 0u; i < transfer_size; ++i)
						this->m_p(i + m_MPISpecData.il, m_MPISpecData.jt + 1) = recv_buf[i];
				});
				buffervector.emplace_back(std::move(send_buffer));
				buffervector.emplace_back(std::move(recv_buffer));
			}

			if(MPI_Waitall(static_cast<int>(used_requests.size()), used_requests.data(), MPI_STATUSES_IGNORE) != MPI_SUCCESS)
				throw std::runtime_error{"MPI_Waitall failed!"};
			for(const auto& resp_func : packaged_responses)
				resp_func();
		}

		void ExchangeVelocityData()
		{
			std::vector<MPI_Request> used_requests;
			std::vector<std::function<void()>> packaged_responses;
			std::vector<std::unique_ptr<T[]>> buffervector;

			if(m_MPISpecData.il != 1u) {
				int rank_source, rank_dest;
				if(MPI_Cart_shift(m_MPIInfo.commhandle, 0, -1, &rank_source, &rank_dest) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Cart_shift failed!"};

				const size_type transfer_size = 2 * ( m_MPISpecData.jt - m_MPISpecData.jb + 1 );
				auto send_buffer = std::make_unique<T[]>(transfer_size);
				auto recv_buffer = std::make_unique<T[]>(transfer_size);
				for(size_type j = 0u; j < transfer_size/2; ++j) {
					send_buffer.get()[j] = m_u(m_MPISpecData.il, j + m_MPISpecData.jb);
					send_buffer.get()[j + transfer_size/2] = m_v(m_MPISpecData.il, j + m_MPISpecData.jb);
				}
				
				MPI_Request send_request, recv_request;
				if(MPI_Isend(send_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::UVRight ), m_MPIInfo.commhandle, &send_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Isend failed!"};
				if(MPI_Irecv(recv_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::UVLeft ), m_MPIInfo.commhandle, &recv_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Irecv failed!"};

				used_requests.emplace_back(std::move(send_request));
				used_requests.emplace_back(std::move(recv_request));
				packaged_responses.emplace_back([recv_buf = recv_buffer.get(), this, transfer_size]() { 
					for(size_type j = 0u; j < transfer_size/2; ++j) {
						this->m_u(m_MPISpecData.il - 1, j + m_MPISpecData.jb) = recv_buf[j];
						this->m_v(m_MPISpecData.il - 1, j + m_MPISpecData.jb) = recv_buf[j + transfer_size/2];
					}
				});
				buffervector.emplace_back(std::move(send_buffer));
				buffervector.emplace_back(std::move(recv_buffer));
			}

			if(m_MPISpecData.ir != m_ProblemGeometry.imax) {
				int rank_source, rank_dest;
				if(MPI_Cart_shift(m_MPIInfo.commhandle, 0, 1, &rank_source, &rank_dest) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Cart_shift failed!"};

				const size_type transfer_size = 2 * ( m_MPISpecData.jt - m_MPISpecData.jb + 1 );
				auto send_buffer = std::make_unique<T[]>(transfer_size);
				auto recv_buffer = std::make_unique<T[]>(transfer_size);
				for(size_type j = 0u; j < transfer_size/2; ++j) {
					send_buffer.get()[j] = m_u(m_MPISpecData.ir, j + m_MPISpecData.jb);
					send_buffer.get()[j + transfer_size/2] = m_v(m_MPISpecData.ir, j + m_MPISpecData.jb);
				}
				
				MPI_Request send_request, recv_request;
				if(MPI_Isend(send_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::UVLeft ), m_MPIInfo.commhandle, &send_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Isend failed!"};
				if(MPI_Irecv(recv_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::UVRight ), m_MPIInfo.commhandle, &recv_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Irecv failed!"};

				used_requests.emplace_back(std::move(send_request));
				used_requests.emplace_back(std::move(recv_request));
				packaged_responses.emplace_back([recv_buf = recv_buffer.get(), this, transfer_size]() { 
					for(size_type j = 0u; j < transfer_size/2; ++j) {
						this->m_u(m_MPISpecData.ir + 1, j + m_MPISpecData.jb) = recv_buf[j];
						this->m_v(m_MPISpecData.ir + 1, j + m_MPISpecData.jb) = recv_buf[j + transfer_size/2];
					}
				});
				buffervector.emplace_back(std::move(send_buffer));
				buffervector.emplace_back(std::move(recv_buffer));
			}

			if(m_MPISpecData.jb != 1u) {
				int rank_source, rank_dest;
				if(MPI_Cart_shift(m_MPIInfo.commhandle, 1, -1, &rank_source, &rank_dest) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Cart_shift failed!"};

				const size_type transfer_size = 2 * ( m_MPISpecData.ir - m_MPISpecData.il + 1 );
				auto send_buffer = std::make_unique<T[]>(transfer_size);
				auto recv_buffer = std::make_unique<T[]>(transfer_size);
				for(size_type i = 0u; i < transfer_size/2; ++i) {
					send_buffer.get()[i] = m_u(i + m_MPISpecData.il, m_MPISpecData.jb);
					send_buffer.get()[i + transfer_size/2] = m_v(i + m_MPISpecData.il, m_MPISpecData.jb);
				}
				
				MPI_Request send_request, recv_request;
				if(MPI_Isend(send_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::UVTop ), m_MPIInfo.commhandle, &send_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Isend failed!"};
				if(MPI_Irecv(recv_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::UVBottom ), m_MPIInfo.commhandle, &recv_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Irecv failed!"};

				used_requests.emplace_back(std::move(send_request));
				used_requests.emplace_back(std::move(recv_request));
				packaged_responses.emplace_back([recv_buf = recv_buffer.get(), this, transfer_size]() { 
					for(size_type i = 0u; i < transfer_size/2; ++i) {
						this->m_u(i + m_MPISpecData.il, m_MPISpecData.jb - 1) = recv_buf[i];
						this->m_v(i + m_MPISpecData.il, m_MPISpecData.jb - 1) = recv_buf[i + transfer_size/2];
					}
				});
				buffervector.emplace_back(std::move(send_buffer));
				buffervector.emplace_back(std::move(recv_buffer));
			}

			if(m_MPISpecData.jt != m_ProblemGeometry.jmax) {
				int rank_source, rank_dest;
				if(MPI_Cart_shift(m_MPIInfo.commhandle, 1, 1, &rank_source, &rank_dest) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Cart_shift failed!"};

				const size_type transfer_size = 2 * ( m_MPISpecData.ir - m_MPISpecData.il + 1 );
				auto send_buffer = std::make_unique<T[]>(transfer_size);
				auto recv_buffer = std::make_unique<T[]>(transfer_size);
				for(size_type i = 0u; i < transfer_size/2; ++i) {
					send_buffer.get()[i] = m_u(i + m_MPISpecData.il, m_MPISpecData.jt);
					send_buffer.get()[i + transfer_size/2] = m_v(i + m_MPISpecData.il, m_MPISpecData.jt);
				}
				
				MPI_Request send_request, recv_request;
				if(MPI_Isend(send_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::UVBottom ), m_MPIInfo.commhandle, &send_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Isend failed!"};
				if(MPI_Irecv(recv_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::UVTop ), m_MPIInfo.commhandle, &recv_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Irecv failed!"};

				used_requests.emplace_back(std::move(send_request));
				used_requests.emplace_back(std::move(recv_request));
				packaged_responses.emplace_back([recv_buf = recv_buffer.get(), this, transfer_size]() { 
					for(size_type i = 0u; i < transfer_size/2; ++i) {
						this->m_u(i + m_MPISpecData.il, m_MPISpecData.jt + 1) = recv_buf[i];
						this->m_v(i + m_MPISpecData.il, m_MPISpecData.jt + 1) = recv_buf[i + transfer_size/2];
					}
				});
				buffervector.emplace_back(std::move(send_buffer));
				buffervector.emplace_back(std::move(recv_buffer));
			}

			if(MPI_Waitall(static_cast<int>(used_requests.size()), used_requests.data(), MPI_STATUSES_IGNORE) != MPI_SUCCESS)
				throw std::runtime_error{"MPI_Waitall failed!"};
			for(const auto& resp_func : packaged_responses)
				resp_func();
		}

		void ExchangeFGData()
		{
			std::vector<MPI_Request> used_requests;
			std::vector<std::function<void()>> packaged_responses;
			std::vector<std::unique_ptr<T[]>> buffervector;

			if(m_MPISpecData.il != 1u) {
				int rank_source, rank_dest;
				if(MPI_Cart_shift(m_MPIInfo.commhandle, 0, -1, &rank_source, &rank_dest) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Cart_shift failed!"};

				const size_type transfer_size = m_MPISpecData.jt - m_MPISpecData.jb + 1;
				auto recv_buffer = std::make_unique<T[]>(transfer_size);
				
				MPI_Request recv_request;
				if(MPI_Irecv(recv_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::F ), m_MPIInfo.commhandle, &recv_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Irecv failed!"};

				used_requests.emplace_back(std::move(recv_request));
				packaged_responses.emplace_back([recv_buf = recv_buffer.get(), this, transfer_size]() { 
					for(size_type j = 0u; j < transfer_size; ++j)
						this->m_F(m_MPISpecData.il - 1, j + m_MPISpecData.jb) = recv_buf[j];
				});
				buffervector.emplace_back(std::move(recv_buffer));
			}

			if(m_MPISpecData.ir != m_ProblemGeometry.imax) {
				int rank_source, rank_dest;
				if(MPI_Cart_shift(m_MPIInfo.commhandle, 0, 1, &rank_source, &rank_dest) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Cart_shift failed!"};

				const size_type transfer_size = m_MPISpecData.jt - m_MPISpecData.jb + 1;
				auto send_buffer = std::make_unique<T[]>(transfer_size);
				for(size_type j = 0u; j < transfer_size; ++j)
					send_buffer.get()[j] = m_F(m_MPISpecData.ir, j + m_MPISpecData.jb);
				
				MPI_Request send_request;
				if(MPI_Isend(send_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::F ), m_MPIInfo.commhandle, &send_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Isend failed!"};

				used_requests.emplace_back(std::move(send_request));
				buffervector.emplace_back(std::move(send_buffer));
			}

			if(m_MPISpecData.jb != 1u) {
				int rank_source, rank_dest;
				if(MPI_Cart_shift(m_MPIInfo.commhandle, 1, -1, &rank_source, &rank_dest) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Cart_shift failed!"};

				const size_type transfer_size = m_MPISpecData.ir - m_MPISpecData.il + 1;
				auto recv_buffer = std::make_unique<T[]>(transfer_size);
				
				MPI_Request recv_request;
				if(MPI_Irecv(recv_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::G ), m_MPIInfo.commhandle, &recv_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Irecv failed!"};

				used_requests.emplace_back(std::move(recv_request));
				packaged_responses.emplace_back([recv_buf = recv_buffer.get(), this, transfer_size]() { 
					for(size_type i = 0u; i < transfer_size; ++i)
						this->m_G(i + m_MPISpecData.il, m_MPISpecData.jb - 1) = recv_buf[i];
				});
				buffervector.emplace_back(std::move(recv_buffer));
			}

			if(m_MPISpecData.jt != m_ProblemGeometry.jmax) {
				int rank_source, rank_dest;
				if(MPI_Cart_shift(m_MPIInfo.commhandle, 1, 1, &rank_source, &rank_dest) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Cart_shift failed!"};

				const size_type transfer_size = m_MPISpecData.ir - m_MPISpecData.il + 1;
				auto send_buffer = std::make_unique<T[]>(transfer_size);
				for(size_type i = 0u; i < transfer_size; ++i)
					send_buffer.get()[i] = m_G(i + m_MPISpecData.il, m_MPISpecData.jt);
				
				MPI_Request send_request;
				if(MPI_Isend(send_buffer.get(), transfer_size, MPIHelper::MPIDatatype_t<T>::value(), rank_dest, static_cast<int>( MPITags_t::G ), m_MPIInfo.commhandle, &send_request) != MPI_SUCCESS)
					throw std::runtime_error{"MPI_Isend failed!"};

				used_requests.emplace_back(std::move(send_request));
				buffervector.emplace_back(std::move(send_buffer));
			}

			if(MPI_Waitall(static_cast<int>(used_requests.size()), used_requests.data(), MPI_STATUSES_IGNORE) != MPI_SUCCESS)
				throw std::runtime_error{"MPI_Waitall failed!"};
			for(const auto& resp_func : packaged_responses)
				resp_func();
		}

		VectorType m_G;
		VectorType m_F;

		problem_data m_ProblemData;
		numeric_parameters m_NumericParameters;
		time_parameters m_TimeParameters;
		problem_geometry m_ProblemGeometry;

		VectorType m_u; ///< Contains the vector \f$ u^{(n)} \f$ at the respective index.
									///< Each time step is represented by a \f$ ( m_ProblemData.imax + 2 ) * ( m_ProblemGeometry.jmax + 2 ) \f$ vector.
		VectorType m_v; ///< Contains the vector \f$ v^{(n)} \f$ at the respective index.
									///< Each time step is represented by a \f$ ( m_ProblemData.imax + 2 ) * ( m_ProblemGeometry.jmax + 2 ) \f$ vector.
		VectorType m_p; 
	};

}

#endif