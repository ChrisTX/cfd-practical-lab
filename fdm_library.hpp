#ifndef FDM_LIBRARY_HPP
#define FDM_LIBRARY_HPP

#include <cassert>
#include <cmath>
#include <type_traits>

namespace FDM {
	template<
		typename TVector2D
	> class FDMVectors2D {
	public:
		using value_type = typename TVector2D::value_type;
		using vector_index_t = typename TVector2D::size_type;
		static_assert(std::is_floating_point<value_type>(), "Can only operate on floating point types!");
		static_assert(Type_Safety::has_compatible_call_operator<TVector2D, vector_index_t, vector_index_t>(), "TVector2D has to be callable as u(i,j)");

		explicit FDMVectors2D(const value_type pγ, const value_type pδx, const value_type pδy, const TVector2D& pu, const TVector2D& pv)
		: γ{pγ}, δx{pδx}, δy{pδy}, u(pu), v(pv) { }

		value_type du2_dx(vector_index_t i, vector_index_t j) const
		{

			value_type res_value = std::pow(u(i,j) + u(i+1,j), 2) - std::pow(u(i-1,j) + u(i,j), 2)
								   + γ * ( std::abs( u(i,j) + u(i+1,j) ) * ( u(i,j) - u(i+1,j) ) - std::abs( u(i-1,j) + u(i,j) ) * ( u(i-1,j) - u(i,j) ) );
			return res_value / (4 * δx);
		}

		value_type duv_dy(vector_index_t i, vector_index_t j) const
		{
			value_type res_value = ( v(i,j) + v(i+1,j) ) * ( u(i,j) + u(i,j+1) ) - ( v(i,j-1) + v(i+1,j-1) ) * ( u(i,j-1) + u(i,j) )
								   + γ * ( std::abs( v(i,j) + v(i+1,j) ) * ( u(i,j) - u(i,j+1) ) - std::abs( v(i,j-1) + v(i+1,j-1) ) * ( u(i,j-1) - u(i,j) ) );
			return res_value / (4 * δy);
		}

		value_type d2u_dx2(vector_index_t i, vector_index_t j) const
		{
			return ( u(i+1,j) - 2 * u(i,j) + u(i-1,j) ) / (std::pow(δx, 2));
		}

		value_type d2u_dy2(vector_index_t i, vector_index_t j) const
		{
			return ( u(i,j+1) - 2 * u(i,j) + u(i,j-1) ) / (std::pow(δy, 2));
		}

		value_type duv_dx(vector_index_t i, vector_index_t j) const
		{
			value_type res_value = ( u(i,j) + u(i,j+1) ) * ( v(i,j) + v(i+1,j) ) - ( u(i-1,j) + u(i-1,j+1) ) * ( v(i-1,j) + v(i,j) )
								   + γ * ( std::abs( u(i,j) + u(i,j+1) ) * ( v(i,j) - v(i+1,j) ) - std::abs( u(i-1,j) + u(i-1,j+1) ) * ( v(i-1,j) - v(i,j) ) );
			return res_value / (4 * δx);
		}

		value_type dv2_dy(vector_index_t i, vector_index_t j) const
		{
			value_type res_value = std::pow(v(i,j) + v(i,j+1), 2) - std::pow(v(i,j-1) + v(i,j), 2)
								   + γ * ( std::abs( v(i,j) + v(i,j+1) ) * ( v(i,j) - v(i,j+1) ) - std::abs( v(i,j-1) + v(i,j) ) * ( v(i,j-1) - v(i,j) ) );
			return res_value / (4 * δy);
		}

		value_type d2v_dx2(vector_index_t i, vector_index_t j) const
		{
			return ( v(i+1,j) - 2 * v(i,j) + v(i-1,j) ) / (std::pow(δx, 2));
		}

		value_type d2v_dy2(vector_index_t i, vector_index_t j) const
		{
			return ( v(i,j+1) - 2 * v(i,j) + v(i,j-1) ) / (std::pow(δy, 2));
		}

	private:
		const TVector2D &u, &v;
		const value_type γ, δx, δy;
	};
}

#endif