#ifndef CFD_TRACER_HPP
#define CFD_TRACER_HPP

#include <cstdlib>
#include <stdexcept>

#include <array>
#include <deque>
#include <functional>
#include <utility>

#include "vtk_library.hpp"
#include "cfd_vtk_library.hpp"


namespace VTK {
	class VTKPointPrinter : public VTKPrinter {
	public:
		template<typename VectorT> 
		explicit VTKPointPrinter(const std::string& FileName, const std::string& vtkDescription, VTKFileFormat_t vtkFormat, const VectorT& PointsList) 
		: VTKPrinter(FileName, vtkDescription, std::move(vtkFormat), 0u, PointsList.size())
		{
			using vector_t = typename VectorT::value_type::value_type;
			if(!PointsList.size())
				throw std::invalid_argument{"Cannot create a file without points!"};

			m_OutputFile << "DATASET POLYDATA\nPOINTS " << PointsList.size() << ' ' << VTKTypeName<typename VectorT::value_type::value_type>::value << '\n';
			for(const auto& cur_point : PointsList) {
				EmitData(cur_point[0u], false);
				EmitData(cur_point[1u], false);
				EmitData(vector_t{0}, true);
			}

			using size_type = typename VectorT::size_type;
			m_OutputFile << "VERTICES " << PointsList.size() << ' ' << 2 * PointsList.size() << '\n';
			for(size_type i = 0u; i < PointsList.size(); ++i) {
				EmitData(1, false);
				EmitData(static_cast<int>(i), true);
			}
		}
	};

	class TimedVTKPointPrinter {
	private:
		const std::string& m_FileNameBase;
		std::size_t m_FileCounter;
		const std::string& m_vtkDescription;
		const VTKFileFormat_t m_vtkFormat;

	public:
		explicit TimedVTKPointPrinter(const std::string& FileNameBase, std::size_t initial_counter, const std::string& vtkDescription, VTKFileFormat_t vtkFormat)
		: m_FileNameBase{FileNameBase}, m_FileCounter{initial_counter}, m_vtkDescription{vtkDescription}, m_vtkFormat{vtkFormat}
		{
			if(FileNameBase.empty())
				throw std::invalid_argument{"Invalid base file name!"};
		}

		template<typename VectorT>
		VTKPrinter PrintNext(const VectorT& PointsList)
		{
			std::ostringstream CurrentFileName(m_FileNameBase, std::ios_base::ate);
			CurrentFileName << m_FileCounter++ << ".vtk";
			return VTKPointPrinter(CurrentFileName.str(), m_vtkDescription, m_vtkFormat, PointsList);
		}
	};
}

namespace CFD {
	template<typename T>
	class ParticleTracerBase {
	protected:
		std::deque<std::array<T, 2u>> m_Particles;
		VTK::TimedVTKPointPrinter m_PointPrinter;

	public:
		using size_type = typename std::deque<std::array<T, 2u>>::size_type;

		ParticleTracerBase(VTK::TimedVTKPointPrinter PointPrinter) : m_PointPrinter(std::move(PointPrinter)) { }

		void InsertParticles(size_type N, T x1, T y1, T x2, T y2)
		{
			assert(N);

			for(size_type i = 0u; i < N; ++i) {
				const T px = x1 + static_cast<T>(i + 1u)/(N + 1u) * ( x2 - x1 );
				const T py = y1 + static_cast<T>(i + 1u)/(N + 1u) * ( y2 - y1 );
				m_Particles.emplace_back(std::array<T, 2u>{px, py});
			}
		}

		void AdvanceParticles(const T δt, const std::function<T(T, T)>& IU, const std::function<T(T, T)>& IV)
		{
			for(auto& cur_particle : m_Particles) {
				const auto u_n = IU(cur_particle[0u], cur_particle[1u]);
				const auto v_n = IV(cur_particle[0u], cur_particle[1u]);

				cur_particle[0u] += δt * u_n;
				cur_particle[1u] += δt * v_n;
			}
		}

		void RemoveParticles(const std::function<bool(T, T)>& is_contained)
		{
			for(size_type i = 0u; i < m_Particles.size(); ++i) {
				if(!is_contained(m_Particles[i][0u], m_Particles[i][1u])) {
					m_Particles.erase(m_Particles.begin() + i);
					--i;
				}
			}
		}
	};

	template<typename T>
	class StreamlineTracer : public ParticleTracerBase<T> {
		using size_type = typename ParticleTracerBase<T>::size_type;
		std::array<T, 4u> m_ParticleLine;
		size_type m_ParticleCount;
		T m_δStreak;
		T m_δInject;

		T m_tStreak = 0;
		T m_tInject = 0;

	public:
		StreamlineTracer(VTK::TimedVTKPointPrinter PointPrinter, T δInject, T δStreak, size_type N, T x1, T y1, T x2, T y2) : ParticleTracerBase<T>(std::move(PointPrinter)), m_δStreak(δStreak), m_δInject(δInject), m_ParticleCount(N), m_ParticleLine{x1, y1, x2, y2} { }

		void Draw(const T δt, const std::function<T(T, T)>& IU, const std::function<T(T, T)>& IV, const std::function<bool(T, T)>& is_contained)
		{
			m_tInject += δt;
			m_tStreak += δt;

			if(m_tInject >= m_δInject) {
				this->InsertParticles(m_ParticleCount, m_ParticleLine[0u], m_ParticleLine[1u], m_ParticleLine[2u], m_ParticleLine[3u]);
				m_tInject = 0;
			}
			if(m_tStreak >= m_δStreak) {
				this->AdvanceParticles(m_tStreak, IU, IV);
				this->RemoveParticles(is_contained);
				this->m_PointPrinter.PrintNext(this->m_Particles);
				m_tStreak = 0;
			}
		}
	};

	template<typename T>
	class ParticleTracer : public ParticleTracerBase<T> {
		using size_type = typename ParticleTracerBase<T>::size_type;
		T m_δTrace;

		T m_tTrace = 0;

	public:
		ParticleTracer(VTK::TimedVTKPointPrinter PointPrinter, T δTrace, size_type N, T x1, T y1, T x2, T y2) : ParticleTracerBase<T>(std::move(PointPrinter)), m_δTrace(δTrace)
		{
			this->InsertParticles(N, x1, y1, x2, y2);
		}

		void Draw(const T δt, const std::function<T(T, T)>& IU, const std::function<T(T, T)>& IV, const std::function<bool(T, T)>& is_contained)
		{
			m_tTrace += δt;

			if(m_tTrace > m_δTrace) {
				this->AdvanceParticles(m_tTrace, IU, IV);
				this->RemoveParticles(is_contained);
				this->m_PointPrinter.PrintNext(this->m_Particles);
				m_tTrace = 0;
			}
		}
	};
}

#endif