#ifndef NUM_VECTOR_HPP
#define NUM_VECTOR_HPP

#include <vector>

namespace numerics {

	namespace {

		template<
			typename T
		> class numerical_vector_operations {
			friend T& operator+=(T&& x, T const& y) { }

			friend T operator+(T&& x, T const& y) { return (x += y); }
			friend T operator+(T const& x, T&& y) { return (x += y); }
			friend T operator+(T const& x, T const& y) { return (T{x} += y); }

			friend T operator-(T&& x, T const& y) { return (x -= y); }
			friend T operator-(T const& x, T&& y) { return (x -= y); }
			friend T operator-(T const& x, T const& y) { return (T{x} -= y); }
		};

	}

	template<
		typename T,
		typename Allocator = std::allocator<T>
	> class numerical_vector final : public std::vector<T, Allocator> {
		using std::vector<T, Allocator>::vector;
	};
}

#endif