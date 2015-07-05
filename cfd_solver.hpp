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

#include "cfd_vtk_library.hpp"
#include "cfd_tracer_library.hpp"
#include "fdm_library.hpp"
#include "linalg_vector2d.hpp"

#include <iostream>

namespace CFD {

	template<
		typename T,
		template<typename...> class VectorT = LinearAlgebra::Vector2D
	> class CFDSolver {
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

		CFDSolver(problem_data pbdata, problem_geometry pbgeom, time_parameters timeparam, numeric_parameters numparam, VectorType u0, VectorType v0, VectorType p0)
		: m_ProblemData(std::move(pbdata)), m_ProblemGeometry(std::move(pbgeom)), m_TimeParameters(std::move(timeparam)), m_NumericParameters(std::move(numparam)), m_u(std::move(u0)), m_v(std::move(v0)), m_p(std::move(p0))
		{
			if(m_u.get_columns() != ( m_ProblemGeometry.imax + 2 ) || m_u.get_rows() != ( m_ProblemGeometry.jmax + 2 ) 
				|| m_u.get_columns() != m_v.get_columns() || m_u.get_rows() != m_v.get_rows() || m_u.get_rows() != m_p.get_rows() || m_u.get_columns() != m_p.get_columns())
				throw std::invalid_argument{"Invalid size for u0, p0 or v0"};
		}

		void Solve()
		{
			T t_since_last_print = 0;
			VTK::TimedVTKRectilinearGridPrinter<T> CFD_VTKPrinters("output-", 0u, "CFD solution", VTK::VTKFileFormat_t::ASCII, static_cast<VTK::VTKIndex_t>(m_ProblemGeometry.imax), static_cast<VTK::VTKIndex_t>(m_ProblemGeometry.jmax), δx(), δy() );
			while(m_TimeParameters.t < m_TimeParameters.t_end) {
				ApplyBoundaryConditions();
				ApplyObstacleConditions();
				SolveLoop();
				t_since_last_print += δt;
				if(t_since_last_print >= m_TimeParameters.t_vis) {
					t_since_last_print = 0;

					auto& vtkprt = CFD_VTKPrinters.GetNextPrinter();
					vtkprt.UsePointData();
					if(!m_ContainedCells.empty()) {
						vtkprt.AddScalars(m_ContainedCells, "Pressure", m_p);
						vtkprt.AddVectors(m_ContainedCells, "Velocity", m_u, m_v);
					}
					else {
						vtkprt.AddScalars("Pressure", m_p);
						vtkprt.AddVectors("Velocity", m_u, m_v);
					}
					// Stream function/vorticity
					VectorType ψ, ζ;
					Computeψ(ψ);
					Computeζ(ζ);
					vtkprt.AddScalars("Stream function", ψ);
					vtkprt.AddScalars("Vorticity", ζ);
				}

				// Streamlines/particles
				auto streamlines_ptr = m_Streamlines.get();
				auto particles_ptr = m_Particles.get();
				if(streamlines_ptr != nullptr) {
					streamlines_ptr->Draw(δt,
						[this](T x, T y) -> auto { return this->InterpolatedU(x, y); },
						[this](T x, T y) -> auto { return this->InterpolatedV(x, y); },
						[this](T x, T y) -> auto { return this->IsContained(x, y); }
					);
				}
				if(particles_ptr != nullptr) {
					particles_ptr->Draw(δt,
						[this](T x, T y) -> auto { return this->InterpolatedU(x, y); },
						[this](T x, T y) -> auto { return this->InterpolatedV(x, y); },
						[this](T x, T y) -> auto { return this->IsContained(x, y); }
					);
				}
			}
		}

		void Computeψ(VectorType& ψ)
		{
			ψ.resize(m_ProblemGeometry.imax + 1, m_ProblemGeometry.jmax + 1);

			for(size_type i = 0u; i <= m_ProblemGeometry.imax; ++i) {
				ψ(i, 0) = 0;
				for(size_type j = 1u; j <= m_ProblemGeometry.jmax; ++j) {
					ψ(i, j) = ψ(i, j - 1);
					if(m_ContainedCells(i, j))
						ψ(i, j) += m_u(i, j) * δx();
				}
			}
		}

		void Computeζ(VectorType& ζ)
		{
			ζ.resize(m_ProblemGeometry.imax + 1, m_ProblemGeometry.jmax + 1, 0);

			for(size_type i = 1u; i < m_ProblemGeometry.imax; ++i)
				for(size_type j = 1u; j < m_ProblemGeometry.jmax; ++j)
					if(m_ContainedCells(i, j) && m_ContainedCells(i + 1, j) && m_ContainedCells(i, j + 1))
						ζ(i, j) = (m_u(i, j + 1) - m_u(i, j)) / δy() - (m_v(i + 1, j) - m_v(i, j)) / δx();
		}

		void SetObstacles(VectorT<bool> contained_cells)
		{
			if(contained_cells.get_rows() != m_ProblemGeometry.jmax + 1 || contained_cells.get_columns() != m_ProblemGeometry.imax + 1)
				throw std::invalid_argument{"Cells mismatch with geometry!"};

			for(size_type i = 1u; i <= m_ProblemGeometry.imax; ++i) {
				for(size_type j = 1u; j <= m_ProblemGeometry.jmax; ++j) {
					if(i == 11 && j == 12)
						printf("I got u now");
					if(!contained_cells(i, j)) {
						using raw_t = std::underlying_type_t<ObstacleCellType_t>;
						raw_t current_cell_type{0};

						if(j && contained_cells(i, j - 1))
							current_cell_type |= static_cast<raw_t>(ObstacleCellType_t::south);
						if(i && contained_cells(i - 1, j))
							current_cell_type |= static_cast<raw_t>(ObstacleCellType_t::west);
						if(j < m_ProblemGeometry.jmax && contained_cells(i, j + 1))
							current_cell_type |= static_cast<raw_t>(ObstacleCellType_t::north);
						if(i < m_ProblemGeometry.imax && contained_cells(i + 1, j))
							current_cell_type |= static_cast<raw_t>(ObstacleCellType_t::east);

						m_ObstacleCells[static_cast<ObstacleCellType_t>(current_cell_type)].emplace(i, j);
					}
				}
			}

			m_ContainedCells = std::move(contained_cells);
			printobs();
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

		T InterpolatedU(T x, T y) const
		{
			const size_type i_star = static_cast<size_type>( x/δx() ) + 1u;
			const size_type j_star = static_cast<size_type>( ( y + δy()/2 )/δy() ) + 1u;

			const T x1 = (i_star - 1) * δx();
			const T x2 = i_star * δx();
			const T y1 = ((j_star - 1) - T{1}/T{2}) * δy();
			const T y2 = (j_star - T{1}/T{2}) * δy();

			return T{1}/( δx() * δy() ) * ( (x2 - x) * (y2 - y) * m_u(i_star - 1, j_star - 1) + (x - x1) * (y2 - y) * m_u(i_star, j_star -1)
				+ (x2 - x) * (y - y1) * m_u(i_star-1, j_star) + (x - x1) * (y - y1) * m_u(i_star, j_star) );
		}

		T InterpolatedV(T x, T y) const
		{
			const size_type i_star = static_cast<size_type>( ( x + δx()/2 )/δx() ) + 1u;
			const size_type j_star = static_cast<size_type>( y/δy() ) + 1u;

			const T x1 = ( (i_star - 1) - T{1}/T{2} ) * δx();
			const T x2 = ( i_star - T{1}/T{2} ) * δx();
			const T y1 = (j_star - 1) * δy();
			const T y2 = j_star * δy();

			return T{1}/( δx() * δy() ) * ( (x2 - x) * (y2 - y) * m_v(i_star - 1, j_star - 1) + (x - x1) * (y2 - y) * m_v(i_star, j_star -1)
				+ (x2 - x) * (y - y1) * m_v(i_star-1, j_star) + (x - x1) * (y - y1) * m_v(i_star, j_star) );
		}

		bool IsContained(T x, T y) const
		{
			const size_type i = static_cast<size_type>( x/δx() );
			const size_type j = static_cast<size_type>( y/δy() );

			if(!i || i > m_ProblemGeometry.imax || !j || j > m_ProblemGeometry.jmax)
				return false;

			auto missing_cells_p = m_ObstacleCells.find(ObstacleCellType_t::notcontained);
			if(missing_cells_p != m_ObstacleCells.end()) {
				auto test_pair = missing_cells_p->second.find(std::make_pair(i, j));
				return (test_pair == missing_cells_p->second.end());
			}

			return true;
		}

		template<class... Args>
		void EnableStreaklines(Args&&... args)
		{
			m_Streamlines = std::make_unique<StreamlineTracer<T>>(std::forward<Args>(args)...);
		}

		template<class... Args>
		void EnableParticles(Args&&... args)
		{
			m_Particles = std::make_unique<ParticleTracer<T>>(std::forward<Args>(args)...);
		}

	private:
		std::unique_ptr<StreamlineTracer<T>> m_Streamlines;
		std::unique_ptr<ParticleTracer<T>> m_Particles;

		T δt; ///< time step size
		inline T δx() const { return m_ProblemGeometry.xlength/m_ProblemGeometry.imax; } 
		inline T δy() const { return m_ProblemGeometry.ylength/m_ProblemGeometry.jmax; }

		void Updateδt()
		{
			if(m_TimeParameters.τ < 0) {
				δt = m_TimeParameters.δt;
				return;
			}

			T umax = 0, vmax = 0;
			for(size_type i = 1u; i <= m_ProblemGeometry.imax; ++i) {
				for(size_type j = 1u; j <= m_ProblemGeometry.jmax; ++j) {
					umax = std::max( umax, std::abs(m_u(i, j)) );
					vmax = std::max( vmax, std::abs(m_v(i, j)) );
				}
			}

			δt = m_TimeParameters.τ * std::min( {
				( m_ProblemData.Re / T{2} ) * std::pow( T{1} / std::pow(δx(), 2) + T{1} / std::pow(δy(), 2), -1),
				δx() / umax,
				δy() / vmax
			} );
		}

		void CalculateFG()
		{
			m_F.resize(m_ProblemGeometry.imax + 1, m_ProblemGeometry.jmax + 1, std::numeric_limits<T>::signaling_NaN());
			m_G.resize(m_ProblemGeometry.imax + 1, m_ProblemGeometry.jmax + 1, std::numeric_limits<T>::signaling_NaN());

			const auto FDM = FDM::FDMVectors2D<VectorType>(m_NumericParameters.γ, δx(), δy(), m_u, m_v);

			for(size_type i = 0u; i <= m_ProblemGeometry.imax; ++i) {
				for(size_type j = 0u; j <= m_ProblemGeometry.jmax; ++j) {
					if(j) {
						if(m_ContainedCells.empty() || !i || i == m_ProblemGeometry.imax || ( m_ContainedCells(i, j) && m_ContainedCells(i + 1, j) ) ) {
							m_F(i, j) = m_u(i, j);

							if(i && i < m_ProblemGeometry.imax)
								m_F(i, j) += δt * ( ( T{1} / m_ProblemData.Re ) * ( FDM.d2u_dx2(i, j) + FDM.d2u_dy2(i, j) ) - FDM.du2_dx(i, j) - FDM.duv_dy(i, j) + m_ProblemData.gx );
							
							assert(std::isfinite(m_F(i, j)));
						}
					}
					if(i) {
						if(m_ContainedCells.empty() || !j || j == m_ProblemGeometry.jmax || ( m_ContainedCells(i, j) && m_ContainedCells(i, j + 1) ) ) {
							m_G(i, j) = m_v(i, j);

							if(j && j < m_ProblemGeometry.jmax)
								m_G(i, j) += δt * ( ( T{1} / m_ProblemData.Re ) * ( FDM.d2v_dx2(i, j) + FDM.d2v_dy2(i, j) ) - FDM.duv_dx(i, j) - FDM.dv2_dy(i, j) + m_ProblemData.gy );
							
							assert(std::isfinite(m_G(i, j)));
						}
					}
				}
			}

			ApplyFGObstacleConditions();
		}

		VectorType CalculateRHS() const
		{
			VectorType RHS(m_ProblemGeometry.imax + 2, m_ProblemGeometry.jmax + 2, std::numeric_limits<T>::signaling_NaN());

			for(size_type i = 1u; i <= m_ProblemGeometry.imax; ++i) {
				for(size_type j = 1u; j <= m_ProblemGeometry.jmax; ++j) {
					if(!m_ContainedCells.empty() && !m_ContainedCells(i, j))
							continue;

					RHS(i, j) = (T{1} / δt) * ( ( m_F(i, j) - m_F(i - 1, j) )/δx() + ( m_G(i, j) - m_G(i, j - 1) )/δy() );
					assert(std::isfinite(RHS(i, j)));
				}
			}

			return RHS;
		}

		void SolveLoop()
		{
			Updateδt();
			CalculateFG();
			const auto RHS = CalculateRHS();

			const auto calc_residual = [&]() -> auto {
				T residual{0};
				for(size_type i = 1u; i <= m_ProblemGeometry.imax; ++i) {
					for(size_type j = 1u; j <= m_ProblemGeometry.jmax; ++j) {
						if(!m_ContainedCells.empty() && !m_ContainedCells(i, j))
							continue;

						residual += std::pow( ( m_p(i+1, j) - 2 * m_p(i, j) + m_p(i-1, j) )/( std::pow(δx(), 2) )
											+ ( m_p(i, j+1) - 2 * m_p(i, j) + m_p(i, j-1) )/( std::pow(δy(), 2) )
											- RHS(i, j), 2);
					}
				}
				assert(std::isfinite(residual));
				return std::sqrt( residual/(m_ProblemGeometry.imax * m_ProblemGeometry.jmax) );
			};

			for(size_type it_count = 0u; it_count < m_NumericParameters.itermax; ++it_count) {
				ApplyPressureConditions();

				for(size_type i = 1u; i <= m_ProblemGeometry.imax; ++i) {
					for(size_type j = 1u; j <= m_ProblemGeometry.jmax; ++j) {
						if(!m_ContainedCells.empty() && !m_ContainedCells(i, j))
							continue;

						m_p(i,j) = (1 - m_NumericParameters.ω) * m_p(i, j) + ( m_NumericParameters.ω / ( T{2} / std::pow(δx(), 2) + T{2} / std::pow(δy(), 2) ) ) *
								 ( ( m_p(i+1,j) + m_p(i-1, j) ) / std::pow(δx(), 2) + ( m_p(i,j+1) + m_p(i, j-1) ) / std::pow(δy(), 2) - RHS(i, j) );
						assert(std::isfinite(m_p(i,j)));
					}
				}

				if(calc_residual() < m_NumericParameters.ε)
					break;
			}

			std::cout << "RES: " << calc_residual() << std::endl;

			for(size_type i = 1u; i < m_ProblemGeometry.imax; ++i)
				for(size_type j = 1u; j <= m_ProblemGeometry.jmax; ++j)
					if(m_ContainedCells.empty() || ( m_ContainedCells(i, j) && m_ContainedCells(i + 1, j) ) )
						m_u(i,j) = m_F(i,j) - δt/δx() * (m_p(i+1,j) - m_p(i,j));

			for(size_type i = 1u; i <= m_ProblemGeometry.imax; ++i)
				for(size_type j = 1u; j < m_ProblemGeometry.jmax; ++j)
					if(m_ContainedCells.empty() || ( m_ContainedCells(i, j) && m_ContainedCells(i, j + 1) ) )
						m_v(i,j) = m_G(i,j) - δt/δy() * (m_p(i,j+1) - m_p(i,j));

			m_TimeParameters.t += δt;
		}

		std::map<BoundaryEdge_t, std::pair<BoundaryCondition_t, input_func_type>> m_BCs;

		enum class ObstacleCellType_t {
			north = 0b0001,
			northeast = 0b0011,
			east = 0b0010,
			southeast = 0b0110,
			south = 0b0100,
			southwest = 0b1100,
			west = 0b1000,
			northwest = 0b1001,
			notcontained = 0b0000,
			interior = 0b1111
		};

		std::map<ObstacleCellType_t, std::set<std::pair<size_type, size_type>>> m_ObstacleCells;
		VectorT<bool> m_ContainedCells;

		void ApplyBoundaryConditions()
		{
			for(std::size_t i = 1; i <= m_ProblemGeometry.imax; ++i) {
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

			for(std::size_t i = 1; i <= m_ProblemGeometry.imax; ++i) {
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

			for(std::size_t j = 1; j <= m_ProblemGeometry.jmax; ++j) {
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

			for(std::size_t j = 1; j <= m_ProblemGeometry.jmax; ++j) {
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

		void ApplyObstacleConditions()
		{
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::north]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_v(i, j) = 0;
				m_u(i - 1, j) = -m_u(i - 1, j + 1);
				m_u(i, j) = -m_u(i, j + 1);
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::south]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_v(i, j - 1) = 0;
				m_u(i - 1, j) = -m_u(i - 1, j - 1);
				m_u(i, j) = -m_u(i, j - 1);
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::west]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_u(i - 1, j) = 0;
				m_v(i, j - 1) = -m_v(i - 1, j - 1);
				m_v(i, j) = -m_v(i - 1, j);
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::east]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_u(i, j) = 0;
				m_v(i, j - 1) = -m_v(i + 1, j - 1);
				m_v(i, j) = -m_v(i + 1, j);
			}

			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::northeast]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_v(i, j) = 0;
				m_u(i, j) = 0;
				m_u(i - 1, j) = -m_u(i - 1, j + 1);
				m_v(i, j - 1) = -m_v(i + 1, j - 1);
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::southeast]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_v(i, j - 1) = 0;
				m_u(i, j) = 0;
				m_u(i - 1, j) = -m_u(i - 1, j - 1);
				m_v(i, j) = -m_v(i + 1, j);
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::southwest]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_v(i, j - 1) = 0;
				m_u(i - 1, j) = 0;
				m_u(i, j) = -m_u(i, j - 1);
				m_v(i, j) = -m_v(i - 1, j);
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::northwest]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_v(i, j) = 0;
				m_u(i - 1, j) = 0;
				m_u(i, j) = -m_u(i, j + 1);
				m_v(i, j - 1) = -m_v(i - 1, j - 1);
			}
		}

		void printobs()
		{
			std::cout << "north:\n";
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::north]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				std::cout << i << ' ' << j << std::endl;
			}
			std::cout << "south:\n";
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::south]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				std::cout << i << ' ' << j << std::endl;
			}
			std::cout << "west:\n";
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::west]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				std::cout << i << ' ' << j << std::endl;
			}
			std::cout << "east:\n";
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::east]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				std::cout << i << ' ' << j << std::endl;
			}
			std::cout << "northeast:\n";
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::northeast]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				std::cout << i << ' ' << j << std::endl;
			}
			std::cout << "southeast:\n";
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::southeast]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				std::cout << i << ' ' << j << std::endl;
			}
			std::cout << "southwest:\n";
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::southwest]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				std::cout << i << ' ' << j << std::endl;
			}
			std::cout << "northwest:\n";
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::northwest]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				std::cout << i << ' ' << j << std::endl;
			}
			std::cout << "outside:\n";
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::notcontained]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				std::cout << i << ' ' << j << std::endl;
			}
		}

		void ApplyFGObstacleConditions()
		{
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::north]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_G(i, j) = m_v(i, j);
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::south]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_G(i, j - 1) = m_v(i, j - 1);
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::west]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_F(i - 1, j) = m_u(i - 1, j);
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::east]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_F(i, j) = m_u(i, j);
			}

			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::northeast]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_F(i, j) = m_u(i, j);
				m_G(i, j) = m_v(i, j);
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::southeast]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_F(i, j) = m_u(i, j);
				m_G(i, j - 1) = m_v(i, j - 1);
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::southwest]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_F(i - 1, j) = m_u(i - 1, j);
				m_G(i, j - 1) = m_v(i, j - 1);
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::northwest]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_F(i - 1, j) = m_u(i - 1, j);
				m_G(i, j) = m_v(i, j);
			}
		}

		void ApplyPressureConditions()
		{
			for(size_type j = 1u; j <= m_ProblemGeometry.jmax; ++j) {
				m_p(0u, j) = m_p(1u, j);
				m_p(m_ProblemGeometry.imax + 1, j) = m_p(m_ProblemGeometry.imax, j);
			}

			for(size_type i = 1u; i <= m_ProblemGeometry.imax; ++i) {
				m_p(i, 0u) = m_p(i, 1u);
				m_p(i, m_ProblemGeometry.jmax + 1) = m_p(i, m_ProblemGeometry.jmax);
			}

			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::north]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_p(i, j) = m_p(i, j + 1);
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::south]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_p(i, j) = m_p(i, j - 1);
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::west]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_p(i, j) = m_p(i - 1, j);
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::east]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_p(i, j) = m_p(i + 1, j);
			}

			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::northeast]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_p(i, j) = ( T{1} / ( std::pow(δx(), 2) + std::pow(δy(), 2) ) ) * ( std::pow(δx(), 2) * m_p(i, j + 1) + std::pow(δy(), 2) * m_p(i + 1, j) );
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::southeast]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_p(i, j) = ( T{1} / ( std::pow(δx(), 2) + std::pow(δy(), 2) ) ) * ( std::pow(δx(), 2) * m_p(i, j - 1) + std::pow(δy(), 2) * m_p(i + 1, j) );
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::southwest]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_p(i, j) = ( T{1} / ( std::pow(δx(), 2) + std::pow(δy(), 2) ) ) * ( std::pow(δx(), 2) * m_p(i, j - 1) + std::pow(δy(), 2) * m_p(i - 1, j) );
			}
			for(const auto& index_pair : m_ObstacleCells[ObstacleCellType_t::northwest]) {
				const auto i = index_pair.first;
				const auto j = index_pair.second;

				m_p(i, j) = ( T{1} / ( std::pow(δx(), 2) + std::pow(δy(), 2) ) ) * ( std::pow(δx(), 2) * m_p(i, j + 1) + std::pow(δy(), 2) * m_p(i - 1, j) );
			}
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