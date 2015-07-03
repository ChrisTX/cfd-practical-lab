#ifndef LINALG_VECTOR2D
#define LINALG_VECTOR2D

#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace LinearAlgebra {
	/*! An \a m × \a n 2D vector with elements of type \a T.
	 *	
	 *	\tparam T The type of the elements in the 2D vector.
	 *	\tparam vector_type Vector type to use as underlying storage.
	 */
	template<
		class T,
		class vector_type = std::vector<T>
	> class Vector2DBase {
	public:
		using value_type = T;
		using size_type = typename vector_type::size_type;
		using difference_type = typename vector_type::difference_type;
		using reference = typename vector_type::reference;
		using const_reference = typename vector_type::const_reference;
		using pointer = typename vector_type::pointer;
		using const_pointer = typename vector_type::const_pointer;
		using iterator = typename vector_type::iterator;
		using const_iterator = typename vector_type::const_iterator;
		using reverse_iterator = typename vector_type::reverse_iterator;
		using const_reverse_iterator = typename vector_type::const_reverse_iterator;

		static_assert(std::is_same<T, typename vector_type::value_type>::value, "Mismatch in value types");

	protected:
		Vector2DBase() = default;

		explicit Vector2DBase(size_type columns, size_type rows) 
			: m_Columns{columns}, m_Rows{rows}, m_UnderlyingVector(columns * rows) { }

		explicit Vector2DBase(size_type columns, size_type rows, vector_type init_vector) 
			: m_Columns{columns}, m_Rows{rows}, m_UnderlyingVector(std::move(init_vector)) { }

		explicit Vector2DBase(size_type columns, size_type rows, value_type value) 
			: m_Columns{columns}, m_Rows{rows}, m_UnderlyingVector(columns * rows, value) { }

	public:
		auto& get_underlying_vector()
		{
			return m_UnderlyingVector;
		}

		decltype(auto) begin()
		{
			return m_UnderlyingVector.begin();
		}

		decltype(auto) begin() const
		{
			return m_UnderlyingVector.begin();
		}

		decltype(auto) cbegin() const
		{
			return m_UnderlyingVector.cbegin();
		}

		decltype(auto) end()
		{
			return m_UnderlyingVector.end();
		}

		decltype(auto) end() const
		{
			return m_UnderlyingVector.end();
		}

		decltype(auto) cend() const
		{
			return m_UnderlyingVector.cend();
		}

		decltype(auto) rbegin()
		{
			return m_UnderlyingVector.rbegin();
		}

		decltype(auto) rbegin() const
		{
			return m_UnderlyingVector.rbegin();
		}

		decltype(auto) crbegin() const
		{
			return m_UnderlyingVector.crbegin();
		}

		decltype(auto) rend()
		{
			return m_UnderlyingVector.rend();
		}

		decltype(auto) rend() const
		{
			return m_UnderlyingVector.rend();
		}

		decltype(auto) crend() const
		{
			return m_UnderlyingVector.crend();
		}

		void swap(Vector2DBase& other)
		{
			std::swap(this->m_UnderlyingVector, other.m_UnderlyingVector);
			std::swap(this->m_Columns, other.m_Columns);
			std::swap(this->m_Rows, other.m_Rows);
		}

		void resize(size_type columns, size_type rows)
		{
			m_Columns = columns;
			m_Rows = rows;
			m_UnderlyingVector.resize(m_Columns * m_Rows);
		}

		void resize(size_type columns, size_type rows, const value_type& value)
		{
			m_Columns = columns;
			m_Rows = rows;
			m_UnderlyingVector.resize(m_Columns * m_Rows, value);
		}

		decltype(auto) size() const
		{
			return m_UnderlyingVector.size();
		}

		decltype(auto) get_rows() const
		{
			return m_Rows;
		}

		decltype(auto) get_columns() const
		{
			return m_Columns;
		}

		decltype(auto) data()
		{
			return m_UnderlyingVector.data();
		}

		decltype(auto) data() const
		{
			return m_UnderlyingVector.data();
		}

		decltype(auto) empty() const
		{
			return m_UnderlyingVector.empty();
		}

	protected:
		size_type m_Rows = 0u;
		size_type m_Columns = 0u;
		vector_type m_UnderlyingVector;
	};

	template<
		class T,
		class vector_type = std::vector<T>
	> class Vector2D : public Vector2DBase<T, vector_type> {
	public:
		using size_type = typename Vector2DBase<T, vector_type>::size_type;
		using value_type = typename Vector2DBase<T, vector_type>::value_type;
		Vector2D() = default;

		explicit Vector2D(size_type columns, size_type rows) 
			: Vector2DBase<T, vector_type>(columns, rows) { }

		explicit Vector2D(size_type columns, size_type rows, vector_type init_vector) 
			: Vector2DBase<T, vector_type>(columns, rows, std::move(init_vector)) { }

		explicit Vector2D(size_type columns, size_type rows, value_type value) 
			: Vector2DBase<T, vector_type>(columns, rows, value) { }

		decltype(auto) at(size_type i, size_type j)
		{
			if(i >= this->m_Columns || j >= this->m_Rows)
				throw std::out_of_range{"Attempted accessing an invalid entry"};

			return this->m_UnderlyingVector.at(i + j * this->m_Columns);
		}

		decltype(auto) at(size_type i, size_type j) const
		{
			if(i >= this->m_Columns || j >= this->m_Rows)
				throw std::out_of_range{"Attempted accessing an invalid entry"};

			return this->m_UnderlyingVector.at(i + j * this->m_Columns);
		}

		decltype(auto) operator()(size_type i, size_type j)
		{
			assert(i < this->m_Columns && j < this->m_Rows);
			return this->m_UnderlyingVector[i + j * this->m_Columns];
		}

		decltype(auto) operator()(size_type i, size_type j) const
		{
			assert(i < this->m_Columns && j < this->m_Rows);
			return this->m_UnderlyingVector[i + j * this->m_Columns];
		}
	};

	template<
		class T,
		class vector_type = std::vector<T>
	> class OffsetVector2D : public Vector2DBase<T, vector_type> {
	public:
		using size_type = typename Vector2DBase<T, vector_type>::size_type;
		using value_type = typename Vector2DBase<T, vector_type>::value_type;
		OffsetVector2D() = default;

		explicit OffsetVector2D(size_type columns, size_type rows, size_type column_off = 0u, size_type row_off = 0u) 
			: Vector2DBase<T, vector_type>(columns, rows), m_RowOffset{row_off}, m_ColumnOffset{column_off} { }

		explicit OffsetVector2D(size_type columns, size_type rows, size_type column_off, size_type row_off, vector_type init_vector) 
			: Vector2DBase<T, vector_type>(columns, rows, std::move(init_vector)), m_RowOffset{row_off}, m_ColumnOffset{column_off} { }

		explicit OffsetVector2D(size_type columns, size_type rows, size_type column_off, size_type row_off, value_type value) 
			: Vector2DBase<T, vector_type>(columns, rows, value), m_RowOffset{row_off}, m_ColumnOffset{column_off} { }

		decltype(auto) at(size_type i, size_type j)
		{
			i -= m_ColumnOffset;
			j -= m_RowOffset;

			if(i >= this->m_Columns || j >= this->m_Rows)
				throw std::out_of_range{"Attempted accessing an invalid entry"};

			return this->m_UnderlyingVector.at(i + j * this->m_Columns);
		}

		decltype(auto) at(size_type i, size_type j) const
		{
			i -= m_ColumnOffset;
			j -= m_RowOffset;

			if(i >= this->m_Columns || j >= this->m_Rows)
				throw std::out_of_range{"Attempted accessing an invalid entry"};

			return this->m_UnderlyingVector.at(i + j * this->m_Columns);
		}

		decltype(auto) operator()(size_type i, size_type j)
		{
			i -= m_ColumnOffset;
			j -= m_RowOffset;

			assert(i < this->m_Columns && j < this->m_Rows);
			return this->m_UnderlyingVector[i + j * this->m_Columns];
		}

		decltype(auto) operator()(size_type i, size_type j) const
		{
			i -= m_ColumnOffset;
			j -= m_RowOffset;
			
			assert(i < this->m_Columns && j < this->m_Rows);
			return this->m_UnderlyingVector[i + j * this->m_Columns];
		}

		void SetOffset(size_type column_off, size_type row_off)
		{
			m_ColumnOffset = column_off;
			m_RowOffset = row_off;
		}
		
	protected:
		size_type m_RowOffset = 0u;
		size_type m_ColumnOffset = 0u;
	};
}

#endif
