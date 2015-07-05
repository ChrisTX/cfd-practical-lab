#ifndef CFD_VTK_LIBRARY_HPP
#define CFD_VTK_LIBRARY_HPP

#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "vtk_library.hpp"

namespace VTK {
	class VTKRectilinearGridPrinter : public VTKPrinter {
	private:
		VTKIndex_t GetIMax() const
		{
			if(m_vtkDatasetAttributes == VTKDatasetAttributes_t::PointData)
				return m_imax;
			else if(m_vtkDatasetAttributes == VTKDatasetAttributes_t::CellData)
				return m_imax - 1;

			throw std::logic_error{"Attempt to write data without dataset attributes!"};
		}

		VTKIndex_t GetJMax() const
		{
			if(m_vtkDatasetAttributes == VTKDatasetAttributes_t::PointData)
				return m_jmax;
			else if(m_vtkDatasetAttributes == VTKDatasetAttributes_t::CellData)
				return m_jmax - 1;

			throw std::logic_error{"Attempt to write data without dataset attributes!"};
		}

		const VTKIndex_t m_imax;
		const VTKIndex_t m_jmax;

	public:
		template<typename T> 
		explicit VTKRectilinearGridPrinter(const std::string& FileName, const std::string& vtkDescription, VTKFileFormat_t vtkFormat, VTKIndex_t imax, VTKIndex_t jmax, T δx, T δy) 
		: VTKPrinter(FileName, vtkDescription, std::move(vtkFormat), (imax - 1) * (jmax - 1), imax * jmax), m_imax{imax}, m_jmax{jmax}
		{
			if(!imax || !jmax)
				throw std::invalid_argument{"Cannot create an empty grid!"};

			m_OutputFile << "DATASET RECTILINEAR_GRID\nDIMENSIONS " << imax << ' ' << jmax << " 1\n";
			
			m_OutputFile << "X_COORDINATES " << imax << ' ' << VTKTypeName<T>::value << '\n';
			for(VTKIndex_t i = 0u; i < imax; ++i)
				EmitData(δx * i, (i + 1 == imax) );

			m_OutputFile << "\nY_COORDINATES " << jmax << ' ' << VTKTypeName<T>::value << '\n';
			for(VTKIndex_t j = 0u; j < jmax; ++j)
				EmitData(δy * j, (j + 1 == jmax) );

			m_OutputFile << "\nZ_COORDINATES 1 " << VTKTypeName<T>::value << "\n0.0";
		}

		template<
			typename TVector2D
		> void AddScalars(const std::string& DataName, TVector2D p, const std::string& LookupTable = "default")
		{
			using vector_index_t = typename TVector2D::size_type;
			static_assert(Type_Safety::has_compatible_call_operator<TVector2D, vector_index_t, vector_index_t>(), "TVector2D has to be callable as u(i,j)");

			if(p.get_columns() != GetIMax() + 2 || p.get_rows() != GetJMax() + 2)
				throw std::invalid_argument{"Unexpected data size!"};
			
			using value_type = typename TVector2D::value_type;
			m_OutputFile << "SCALARS " << DataName << ' ' << VTKTypeName<value_type>::value << '\n';
			m_OutputFile << "LOOKUP_TABLE " << LookupTable << '\n';

			for(vector_index_t j = 1u; j < GetJMax() + 1u; ++j)
				for(vector_index_t i = 1u; i < GetIMax() + 1u; ++i)
					EmitData(p(i,j), true);
		}

		template<
			typename TVector2D
		> void AddVectors(const std::string& DataName, TVector2D u, TVector2D v)
		{
			using vector_index_t = typename TVector2D::size_type;
			static_assert(Type_Safety::has_compatible_call_operator<TVector2D, vector_index_t, vector_index_t>(), "TVector2D has to be callable as u(i,j)");

			if(u.get_columns() != v.get_columns() || u.get_rows() != v.get_rows() || u.get_columns() != GetIMax() + 2 || u.get_rows() != GetJMax() + 2)
				throw std::invalid_argument{"Unexpected data size!"};
			
			using value_type = typename TVector2D::value_type;
			m_OutputFile << "VECTORS " << DataName << ' ' << VTKTypeName<value_type>::value << '\n';

			for(vector_index_t j = 1u; j < GetJMax() + 1u; ++j) {
				for(vector_index_t i = 1u; i < GetIMax() + 1u; ++i) {
					EmitData(u(i,j), false);
					EmitData(v(i,j), false);
					EmitData(0, true);
				}
			}
		}

		template<
			typename TVector2D,
			typename TTestFunc
		> void AddScalars(const TTestFunc& is_contained, const std::string& DataName, TVector2D p, const std::string& LookupTable = "default")
		{
			using vector_index_t = typename TVector2D::size_type;
			static_assert(Type_Safety::has_compatible_call_operator<TVector2D, vector_index_t, vector_index_t>(), "TVector2D has to be callable as u(i,j)");
			static_assert(Type_Safety::has_compatible_call_operator<TTestFunc, vector_index_t, vector_index_t>()
				&& std::is_same<bool, std::result_of_t<TTestFunc(vector_index_t, vector_index_t)>>(), "TTestFunc has to be callable as F(i,j) and return bool");

			if(p.get_columns() != GetIMax() + 2 || p.get_rows() != GetJMax() + 2)
				throw std::invalid_argument{"Unexpected data size!"};
			
			using value_type = typename TVector2D::value_type;
			m_OutputFile << "SCALARS " << DataName << ' ' << VTKTypeName<value_type>::value << '\n';
			m_OutputFile << "LOOKUP_TABLE " << LookupTable << '\n';

			for(vector_index_t j = 1u; j < GetJMax() + 1u; ++j)
				for(vector_index_t i = 1u; i < GetIMax() + 1u; ++i)
					EmitData(p(i,j), true);
		}

		template<
			typename TVector2D,
			typename TTestFunc
		> void AddVectors(const TTestFunc& is_contained, const std::string& DataName, TVector2D u, TVector2D v)
		{
			using vector_index_t = typename TVector2D::size_type;
			using value_t = typename TVector2D::value_type;
			static_assert(Type_Safety::has_compatible_call_operator<TVector2D, vector_index_t, vector_index_t>(), "TVector2D has to be callable as u(i,j)");
			static_assert(Type_Safety::has_compatible_call_operator<TTestFunc, vector_index_t, vector_index_t>()
				&& std::is_same<bool, std::result_of_t<TTestFunc(vector_index_t, vector_index_t)>>(), "TTestFunc has to be callable as F(i,j) and return bool");

			if(u.get_columns() != v.get_columns() || u.get_rows() != v.get_rows() || u.get_columns() != GetIMax() + 2 || u.get_rows() != GetJMax() + 2)
				throw std::invalid_argument{"Unexpected data size!"};
			
			using value_type = typename TVector2D::value_type;
			m_OutputFile << "VECTORS " << DataName << ' ' << VTKTypeName<value_type>::value << '\n';

			for(vector_index_t j = 1u; j < GetJMax() + 1u; ++j) {
				for(vector_index_t i = 1u; i < GetIMax() + 1u; ++i) {
					if(is_contained(i, j)) {
						EmitData(u(i,j), false);
						EmitData(v(i,j), false);
					}
					else {
						EmitData(value_t{0}, false);
						EmitData(value_t{0}, false);
					}
					EmitData(value_t{0}, true);
				}
			}
		}
	};

	template<
		typename T
	> class TimedVTKRectilinearGridPrinter {
	private:
		const std::string& m_FileNameBase;
		std::size_t m_FileCounter;
		const std::string& m_vtkDescription;
		const VTKFileFormat_t m_vtkFormat;
		const VTKIndex_t m_imax;
		const VTKIndex_t m_jmax;
		T m_δx;
		T m_δy;

		std::unique_ptr<VTKRectilinearGridPrinter> m_curPrinter;

	public:
		explicit TimedVTKRectilinearGridPrinter(const std::string& FileNameBase, std::size_t initial_counter, const std::string& vtkDescription, VTKFileFormat_t vtkFormat, VTKIndex_t imax, VTKIndex_t jmax, T δx, T δy)
		: m_FileNameBase{FileNameBase}, m_FileCounter{initial_counter}, m_vtkDescription{vtkDescription}, m_vtkFormat{vtkFormat}, m_imax{imax}, m_jmax{jmax}, m_δx{δx}, m_δy{δy}
		{
			if(FileNameBase.empty())
				throw std::invalid_argument{"Invalid base file name!"};
		}

		VTKRectilinearGridPrinter& GetPrinter()
		{
			auto pprint = m_curPrinter.get();
			if(!pprint)
				throw std::logic_error{"No printer created!"};
			return *pprint;
		}

		void GenerateNextPrinter()
		{
			std::ostringstream CurrentFileName(m_FileNameBase, std::ios_base::ate);
			CurrentFileName << m_FileCounter++ << ".vtk";
			m_curPrinter = std::make_unique<VTKRectilinearGridPrinter>(CurrentFileName.str(), m_vtkDescription, m_vtkFormat, m_imax, m_jmax, m_δx, m_δy);
		}

		VTKRectilinearGridPrinter& GetNextPrinter()
		{
			GenerateNextPrinter();
			return GetPrinter();
		}
	};
}

#endif
