/*
 * IO.hpp
 *
 *  Created on: Feb 5, 2026
 *      Author: dmarce1
 */

#ifndef IO_HPP_
#define IO_HPP_

#include <cstdlib>
#include <string>

template <typename... Args>
std::string print2string(const char *format, Args... args) {
	static thread_local char *buffer = NULL;
	static thread_local size_t blen = 1;
	size_t const len = snprintf(NULL, 0, format, args...) + 1;
	if (len > blen) {
		if (buffer) {
			free(buffer);
		}
		while (len > blen) {
			blen *= 2;
		}
		buffer = (char *)malloc(sizeof(char) * blen);
	}
	snprintf(buffer, blen, format, args...);
	return std::string(buffer);
}

template <std::integral T>
std::string toBinary(T i) {
	constexpr int S = std::numeric_limits<T>::is_signed;
	constexpr int N = std::numeric_limits<T>::digits + S;
	static thread_local char str[N + 1];
	int k = N;
	str[k--] = '\0';
	if constexpr (S) {
		if (i >= 0) {
			str[0] = '+';
		} else {
			str[0] = '-';
			i = -i;
		}
	}
	while (k >= S) {
		str[k--] = '0' + (i & 1);
		i >>= 1;
	}
	auto const rc = std::string(str);
	return rc;
}

#endif /* IO_HPP_ */
