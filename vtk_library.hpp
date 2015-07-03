#ifndef VTK_LIBRARY_HPP
#define VTK_LIBRARY_HPP

#include <cassert>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "endian_converter.hpp"
#include "type_safety.hpp"

namespace VTK {
	// These dataType names are defined in the VTK file format
	template<
		typename T
	> struct VTKTypeName { };

	template<> struct VTKTypeName<bool> {
		static constexpr auto value = "bit";
	};

	template<> struct VTKTypeName<unsigned char> {
		static constexpr auto value = "unsigned_char";
	};

	template<> struct VTKTypeName<char> {
		static constexpr auto value = "char";
	};

	template<> struct VTKTypeName<unsigned short> {
		static constexpr auto value = "unsigned_short";
	};

	template<> struct VTKTypeName<short> {
		static constexpr auto value = "short";
	};

	template<> struct VTKTypeName<unsigned int> {
		static constexpr auto value = "unsigned_int";
	};

	template<> struct VTKTypeName<int> {
		static constexpr auto value = "int";
	};

	template<> struct VTKTypeName<unsigned long> {
		static constexpr auto value = "unsigned_long";
	};

	template<> struct VTKTypeName<long> {
		static constexpr auto value = "long";
	};

	template<> struct VTKTypeName<float> {
		static constexpr auto value = "float";
	};

	template<> struct VTKTypeName<double> {
		static constexpr auto value = "double";
	};

	enum class VTKFileFormat_t {
		ASCII,
		BINARY
	};

	namespace detail {
		template<
			typename T,
			typename K = decltype(VTKTypeName<T>::value)  
		> constexpr bool act_is_vtk_type(int) { return true; }

		template<
			typename T
		> constexpr bool act_is_vtk_type(float) { return false; }
	}

	template<
		typename T
	> using is_VTK_supported_data_type = std::integral_constant<bool, detail::act_is_vtk_type<T>(0)>;

	// Quote from the VTK documentation "Cell types and indices are of type int."
	using VTKIndex_t = int;

	enum class VTKDatasetAttributes_t {
		Unspecified,
		PointData,
		CellData
	};

	class VTKPrinter {
	public:
		void UsePointData()
		{
			if(!m_Points)
				throw std::logic_error{"No points defined"};

			m_vtkDatasetAttributes = VTKDatasetAttributes_t::PointData;
			m_OutputFile << "POINT_DATA " << m_Points << '\n';
		}

		void UseCellData()
		{
			if(!m_Cells)
				throw std::logic_error{"No cells defined"};

			m_vtkDatasetAttributes = VTKDatasetAttributes_t::CellData;
			m_OutputFile << "CELL_DATA " << m_Cells << '\n';
		}

		VTKDatasetAttributes_t GetDatasetAttributes() const
		{
			return m_vtkDatasetAttributes;
		}

		VTKIndex_t GetNumberOfCells() const
		{
			return m_Cells;
		}

		VTKIndex_t GetNumberOfPoints() const
		{
			return m_Points;
		}

		VTKFileFormat_t GetFileFormat() const
		{
			return m_vtkFormat;
		}

		template<
			typename InputIt
		> void AddScalars(const std::string& DataName, InputIt first, InputIt last, const std::string& LookupTable = "default")
		{
			static_assert(std::is_base_of<std::input_iterator_tag, typename std::iterator_traits<InputIt>::iterator_category>(), "InputIterator required");

			VTKIndex_t DataCounter = GetActiveDataSize();
			
			using value_type = typename std::iterator_traits<InputIt>::value_type;
			m_OutputFile << "SCALARS " << DataName << ' ' << VTKTypeName<value_type>::value << '\n';
			m_OutputFile << "LOOKUP_TABLE " << LookupTable << '\n';
			

			for(auto it = first; it != last; ++it) {
				--DataCounter;
				EmitData(*it, true);
			}

			if(DataCounter != 0)
				throw std::invalid_argument{"Unexpected data size!"};
		}

		template<
			typename InputIt
		> void AddVectors(const std::string& DataName, InputIt x_first, InputIt x_last, InputIt y_first, InputIt y_last, InputIt z_first, InputIt z_last)
		{
			static_assert(std::is_base_of<std::input_iterator_tag, typename std::iterator_traits<InputIt>::iterator_category>(), "InputIterator required");

			VTKIndex_t DataCounter = GetActiveDataSize();
			
			using value_type = typename std::iterator_traits<InputIt>::value_type;
			m_OutputFile << "VECTORS " << DataName << ' ' << VTKTypeName<value_type>::value << '\n';

			for(; x_first != x_last && y_first != y_last && z_first != z_last;) {
				--DataCounter;
				EmitData(*x_first, false);
				EmitData(*y_first, false);
				EmitData(*z_first, true);

				++x_first;
				++y_first;
				++z_first;
			}

			if(DataCounter != 0)
				throw std::invalid_argument{"Unexpected data size!"};
		}

	protected:
		explicit VTKPrinter(const std::string& FileName, const std::string& vtkDescription, VTKFileFormat_t vtkFormat, VTKIndex_t Cells, VTKIndex_t Points) :
		m_Cells{Cells}, m_Points{Points}, m_vtkFormat{vtkFormat}
		{
			// The value of 256 is given in the VTK file format documentation.
			if(vtkDescription.size() > 256u)
				throw std::invalid_argument{"Overlong VTK description: There is a 256 character limit."};

			m_OutputFile.exceptions(std::ios::failbit);
			m_OutputFile.open(FileName, std::ios::binary);
			m_OutputFile << "# vtk DataFile Version 3.0\n" << vtkDescription << '\n';
			switch(vtkFormat) {
				case VTKFileFormat_t::ASCII:
					m_OutputFile << "ASCII\n";
					break;

				case VTKFileFormat_t::BINARY:
					m_OutputFile << "BINARY\n";
					break;
			}
		}

		VTKIndex_t GetActiveDataSize() const
		{
			if(m_vtkDatasetAttributes == VTKDatasetAttributes_t::PointData)
				return m_Points;
			else if(m_vtkDatasetAttributes == VTKDatasetAttributes_t::CellData)
				return m_Cells;

			throw std::logic_error{"Attempt to write data without dataset attributes!"};
		}

		template<
			typename T
		> void EmitData(T value, bool last_in_line, bool force_endline = false)
		{
			static_assert(is_VTK_supported_data_type<T>(), "T is not supported by VTK");

			if(m_vtkFormat == VTKFileFormat_t::ASCII) {
				m_OutputFile << value << (last_in_line ? '\n' : ' ');
			}
			else {
				unsigned char raw_data[sizeof(T)];
				Utility::to_big_endian(value, raw_data);
				m_OutputFile.write(reinterpret_cast<char*>( raw_data ), sizeof(T));
				if(force_endline && last_in_line)
					m_OutputFile << '\n';
			}
		}

		VTKFileFormat_t m_vtkFormat;
		VTKIndex_t m_Cells = 0;
		VTKIndex_t m_Points = 0;
		VTKDatasetAttributes_t m_vtkDatasetAttributes = VTKDatasetAttributes_t::Unspecified;
		std::ofstream m_OutputFile;
	};
}

#endif