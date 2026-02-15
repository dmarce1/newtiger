/*
 * SubArray.hpp
 *
 *  Created on: Feb 13, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_SUBARRAY_HPP_
#define INCLUDE_SUBARRAY_HPP_

#include <array>
#include <cassert>
#include <concepts>
#include <cstdint>
#include <unordered_map>
#include <vector>

template <class T>
concept IntegralArray = requires {
	typename T::value_type;
	std::tuple_size<T>::value;
} && std::integral<typename T::value_type> && std::same_as<T, std::array<typename T::value_type, std::tuple_size_v<T>>>;

template <size_t D>
class Box {
	static constexpr size_t childCount = 1 << D;
	using UInt = std::uint64_t;
	using SInt = std::int64_t;
	std::array<UInt, D> b_{};
	std::array<UInt, D> e_{};

public:
	Box() = default;
	Box(Box const &) = default;
	Box(Box &&) = default;
	Box &operator=(Box const &) = default;
	Box &operator=(Box &&) = default;
	Box(IntegralArray auto const &b, IntegralArray auto const &e) :
		b_{b}, e_{e} {
	}
	Box(UInt n) :
		e_{n} {
	}
	auto const &begin() const {
		return b_;
	}
	auto const &end() const {
		return e_;
	}
	bool contains(IntegralArray auto const &point) const {
		for (size_t d = 0; d < D; ++d) {
			if (point[d] < b_[d]) return false;
			if (point[d] >= e_[d]) return false;
		}
		return true;
	}
	bool contains(Box const &other) const {
		if (!contains(other.b_)) return false;
		if (!contains(other.e_)) return false;
		return true;
	}
	bool intersects(Box const &A) const {
		for (size_t d = 0; d < D; ++d) {
			if (b_[d] >= A.e_[d]) return false;
			if (e_[d] <= A.b_[d]) return false;
		}
		return true;
	}
	bool operator==(Box const &other) const {
		if (b_ != other.b_) return false;
		if (e_ != other.e_) return false;
		return true;
	}
	bool operator!=(Box const &other) const {
		return !(*this == other);
	}
	UInt span(size_t d) const {
		return std::max(e_[d] - b_[d], 0);
	}
	UInt span() const {
		std::array<UInt, D> s;
		for (size_t d = 0; d < D; ++d) {
			s[d] = span(s);
		}
		return s;
	}
	UInt volume() const {
		UInt v = 1;
		for (size_t d = 0; d < D; ++d) {
			v *= span(d);
		}
		return v;
	}
	Box expand(SInt n) const {
		Box b;
		for (size_t d = 0; d < D; ++d) {
			b.b_[d] -= n;
			b.e_[d] += n;
		}
		return b;
	}
	Box shrink(SInt n) const {
		return expand(-n);
	}
	Box scale(UInt n) const {
		Box b;
		for (size_t d = 0; d < D; ++d) {
			b.b_[d] *= n;
			b.e_[d] *= n;
		}
		return b;
	}
	Box shift(IntegralArray auto const &n) const {
		Box b;
		for (size_t d = 0; d < D; ++d) {
			b.b_[d] = b_[d] + n[d];
			b.e_[d] = e_[d] + n[d];
		}
		return b;
	}
	Box shift(SInt n, size_t d) const {
		std::array<SInt, D> s;
		for (size_t i = 0; i < D; ++i) {
			s[i] = (i == d) ? n : 0;
		}
		return shift(s);
	}
	size_t flatten(IntegralArray auto const &I) const {
		size_t i = 0;
		for (int d = 0; d < D; d++) {
			i = span(d) * i + (I[d] - b_[d]);
		}
		return i;
	}
	friend Box bounding(Box const &A, Box const &B) {
		Box I;
		for (size_t d = 0; d < D; ++d) {
			I.begin_[d] = std::min(A.begin_[d], B.begin_[d]);
			I.end_[d] = std::max(A.end_[d], B.end_[d]);
		}
		return I;
	}
	friend Box intersection(Box const &A, Box const &B) {
		Box I;
		for (size_t d = 0; d < D; ++d) {
			I.begin_[d] = std::max(A.begin_[d], B.begin_[d]);
			I.end_[d] = std::min(A.end_[d], B.end_[d]);
		}
		return I;
	}
	void serialize(auto &arc, unsigned) {
		arc & b_;
		arc & e_;
	}
};

template <typename T, size_t D>
class SubArray {
	Box<D> box_;
	std::vector<T> data_;

public:
	using Index = std::array<size_t, D>;
	SubArray() = default;
	SubArray(size_t n) :
		box_{n}, data_{box_.volume()} {
	}
	SubArray(SubArray const &) = default;
	SubArray(SubArray &&) = default;
	SubArray &operator=(SubArray const &) = default;
	SubArray &operator=(SubArray &&) = default;
	SubArray get(Box<D> sbox) const {
		assert(box_.contains(sbox));
		SubArray sub(sbox);
		auto const begin = sbox.begin();
		auto const end = sbox.end();
		auto const count = sbox.volume();
		auto I = begin;
		for (size_t i = 0; i < count; i++) {
			auto const j = box_.flatten(I);
			sub.data_[i] = data_[j];
			for (int d = D - 1; d >= 0; --d) {
				if (++I[d] < end[d]) break;
				I[d] = begin[d];
			}
		}
	}
	SubArray &set(SubArray const &sub) {
		auto const sbox = sub.box_;
		assert(box_.contains(sbox));
		auto const begin = sbox.begin();
		auto const end = sbox.end();
		auto const count = sbox.volume();
		auto I = begin;
		for (size_t i = 0; i < count; i++) {
			auto const j = box_.flatten(I);
			data_[j] = sub.data_[i];
			for (int d = D - 1; d >= 0; --d) {
				if (++I[d] < end[d]) break;
				I[d] = begin[d];
			}
		}
		return *this;
	}
	void serialize(auto &arc, unsigned) {
		arc & box_;
		arc & data_;
	}
};

#endif /* INCLUDE_SUBARRAY_HPP_ */
