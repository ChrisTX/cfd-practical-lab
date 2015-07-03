#ifndef TYPE_SAFETY_HPP
#define TYPE_SAFETY_HPP

#include <type_traits>

namespace Type_Safety {
	inline namespace {
		namespace detail {
			template<
				typename T,
				typename... Args,
				typename K = std::result_of_t<T(Args...)>
			> constexpr bool act_has_comp_call_op(int) { return true; }

			template<
				typename T,
				typename... Args
			> constexpr bool act_has_comp_call_op(float) { return false; }
		}

		template<
			typename T,
			typename... Args
		> using has_compatible_call_operator = std::integral_constant<bool, detail::act_has_comp_call_op<T, Args...>(0)>;

		namespace detail {
			template<
				typename T,
				typename... Args,
				typename K = std::result_of_t<decltype(&T::at)(T, Args...)>
			> constexpr bool act_has_comp_at_mem(int) { return true; }

			template<
				typename T,
				typename... Args
			> constexpr bool act_has_comp_at_mem(float) { return false; }
		}

		template<
			typename T,
			typename... Args
		> using has_compatible_at_member = std::integral_constant<bool, detail::act_has_comp_at_mem<T, Args...>(0)>;
	}
}

#endif