#ifndef ENDIAN_CONVERTER_HPP
#define ENDIAN_CONVERTER_HPP

#include <algorithm>
#include <cstdint>
#include <type_traits>

namespace Utility {
	enum class endianness_t {
		unknown,
		big_endian,
		little_endian_byte,
		little_endian_word,
		PDP_endian
	};

	inline endianness_t determine_endianness() {
		const std::uint32_t ref_value{0x0A0B0C0D};
		const unsigned char* obj_repres = reinterpret_cast<const unsigned char*>(&ref_value);

		if(obj_repres[0] == '\x0A' && obj_repres[1] == '\x0B' && obj_repres[2] == '\x0C' && obj_repres[3] == '\x0D')
			return endianness_t::big_endian;
		else if(obj_repres[0] == '\x0D' && obj_repres[1] == '\x0C' && obj_repres[2] == '\x0B' && obj_repres[3] == '\x0A')
			return endianness_t::little_endian_byte;
		else if(obj_repres[0] == '\x0C' && obj_repres[1] == '\x0D' && obj_repres[2] == '\x0A' && obj_repres[3] == '\x0B')
			return endianness_t::little_endian_word;
		else if(obj_repres[0] == '\x0B' && obj_repres[1] == '\x0A' && obj_repres[2] == '\x0D' && obj_repres[3] == '\x0C')
			return endianness_t::PDP_endian;
		else
			return endianness_t::unknown;
	}

	template<
		typename T
	> inline std::enable_if_t<(sizeof(T) > 1u), void> to_big_endian(const T org_value, unsigned char conv_obj_repres[sizeof(T)]) {
		// is_arithmetic implies TriviallyCopyable, which we need
		static_assert(std::is_arithmetic<T>(), "Only arithmetic types can be converted");
		static_assert(sizeof(T) % 2 == 0u, "Only types with an even size are supported");

		const endianness_t cur_endianness = determine_endianness();
		if(cur_endianness == endianness_t::big_endian) {
			std::copy(reinterpret_cast<const unsigned char*>(&org_value), reinterpret_cast<const unsigned char*>(&org_value) + sizeof(T), conv_obj_repres);
		}
		else if(cur_endianness == endianness_t::little_endian_byte || cur_endianness == endianness_t::little_endian_word) {
			const std::size_t endian_inc = (cur_endianness == endianness_t::little_endian_word) ? 2u : 1u;
			const unsigned char* obj_repres = reinterpret_cast<const unsigned char*>(&org_value);

			for(std::size_t i = 0; i < sizeof(T); i += endian_inc)
				std::copy((obj_repres + sizeof(T)) - (i + 1) * endian_inc, (obj_repres + sizeof(T)) - i * endian_inc, conv_obj_repres + i * endian_inc);
		}
		else {
			throw std::logic_error{"The current endianness is not being supported"};
		}
	}

	template<
		typename T
	> inline std::enable_if_t<sizeof(T) == 1u, void> to_big_endian(const T org_value, unsigned char conv_obj_repres[sizeof(T)]) {
		std::copy(reinterpret_cast<const unsigned char*>(&org_value), reinterpret_cast<const unsigned char*>(&org_value) + sizeof(T), conv_obj_repres);
	}
}

#endif