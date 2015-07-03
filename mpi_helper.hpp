#ifndef MPI_HELPER_HPP
#define MPI_HELPER_HPP

#include <mpi.h>

namespace MPIHelper {
	template<
		typename T
	> struct MPIDatatype_t { };

	template<> struct MPIDatatype_t<char> {
		static decltype(auto) value()
		{
			return MPI_CHAR;
		}
	};

	template<> struct MPIDatatype_t<double> {
		static decltype(auto) value()
		{
			return MPI_DOUBLE;
		}
	};

	template<> struct MPIDatatype_t<float> {
		static decltype(auto) value()
		{
			return MPI_FLOAT;
		}
	};

	template<> struct MPIDatatype_t<int> {
		static decltype(auto) value()
		{
			return MPI_INT;
		}
	};

	template<> struct MPIDatatype_t<long> {
		static decltype(auto) value()
		{
			return MPI_LONG;
		}
	};

	template<> struct MPIDatatype_t<long long> {
		static decltype(auto) value()
		{
			return MPI_LONG_LONG_INT;
		}
	};

	template<> struct MPIDatatype_t<long double> {
		static decltype(auto) value()
		{
			return MPI_LONG_DOUBLE;
		}
	};

	template<> struct MPIDatatype_t<short> {
		static decltype(auto) value()
		{
			return MPI_SHORT;
		}
	};

	template<> struct MPIDatatype_t<unsigned char> {
		static decltype(auto) value()
		{
			return MPI_UNSIGNED_CHAR;
		}
	};

	template<> struct MPIDatatype_t<unsigned int> {
		static decltype(auto) value()
		{
			return MPI_UNSIGNED;
		}
	};

	template<> struct MPIDatatype_t<unsigned long> {
		static decltype(auto) value()
		{
			return MPI_UNSIGNED_LONG;
		}
	};

	template<> struct MPIDatatype_t<unsigned short> {
		static decltype(auto) value()
		{
			return MPI_UNSIGNED_SHORT;
		}
	};
}

#endif