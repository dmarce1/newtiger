/*
 * SubArray.hpp
 *
 *  Created on: Feb 13, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_BOX_HPP_
#define INCLUDE_BOX_HPP_

#include "Integer.hpp"

#include <algorithm>
#include <array>
#include <ostream>

template <Integer D>
class Box {
	static constexpr Integer childCount = 1 << D;
	std::array<Integer, D> begin_;
	std::array<Integer, D> end_;

public:
	constexpr Box() = default;
	constexpr Box(Box const &) = default;
	constexpr Box(Box &&) = default;
	constexpr Box(std::array<Integer, D> const &b, std::array<Integer, D> const &e) {
		std::copy_n(b.begin(), D, begin_.begin());
		std::copy_n(e.begin(), D, end_.begin());
	}
	constexpr Box(Integer e) {
		begin_.fill(0_I);
		end_.fill(e);
	}
	constexpr Box(Integer b, Integer e) {
		begin_.fill(b);
		end_.fill(e);
	}
	constexpr Box &operator=(Box const &) = default;
	constexpr Box &operator=(Box &&) = default;
	constexpr auto const &begin() const {
		return begin_;
	}
	constexpr auto const &end() const {
		return end_;
	}
	constexpr auto span() const {
		std::array<Integer, D> s;
		for (Integer d = 0; d < D; ++d) {
			s[d] = span(d);
		}
		return s;
	}
	constexpr bool contains(std::array<Integer, D> const &point) const {
		for (Integer d = 0; d < D; ++d) {
			if (point[d] < begin_[d]) return false;
			if (point[d] >= end_[d]) return false;
		}
		return true;
	}
	constexpr bool contains(Box const &other) const {
		if (!contains(other.begin_)) return false;
		if (!contains(other.end_)) return false;
		return true;
	}
	constexpr bool intersects(Box const &A) const {
		for (Integer d = 0; d < D; ++d) {
			if (begin_[d] >= A.end_[d]) return false;
			if (end_[d] <= A.begin_[d]) return false;
		}
		return true;
	}
	constexpr bool operator==(Box const &other) const {
		if (begin_ != other.begin_) return false;
		if (end_ != other.end_) return false;
		return true;
	}
	constexpr bool operator!=(Box const &other) const {
		return !(*this == other);
	}
	constexpr bool operator<(Box const &other) const {
		return !(*this == other);
	}
	constexpr Integer span(Integer d) const {
		return std::max(end_[d] - begin_[d], 0_I);
	}
	constexpr Integer volume() const {
		Integer v = 1_I;
		for (Integer d = 0; d < D; ++d) {
			v *= span(d);
		}
		return v;
	}
	constexpr Box expand(Integer n) const {
		Box b;
		for (Integer d = 0; d < D; ++d) {
			b.begin_[d] = begin_[d] - n;
			b.end_[d] = end_[d] + n;
		}
		return b;
	}
	constexpr Box scale(Integer n) const {
		Box b;
		for (Integer d = 0; d < D; ++d) {
			b.begin_[d] = begin_[d] * n;
			b.end_[d] = end_[d] * n;
		}
		return b;
	}
	constexpr Box<D - 1> slice(Integer d) const {
		Box<D - 1> b;
		for (Integer d1 = 0; d1 < D - 1; ++d1) {
			Integer const d2 = (d1 >= d) ? (d1 + 1) : d1;
			b.begin_[d1] = begin_[d2];
			b.end_[d1] = end_[d2];
		}
		return b;
	}
	constexpr Box<D> slice(Integer d, Integer n) const {
		auto b = *this;
		b.begin_[d] = n;
		b.end_[d] = n + 1;
		return b;
	}
	constexpr Box shift(std::array<Integer, D> const &n) const {
		Box b;
		for (Integer d = 0; d < D; ++d) {
			b.begin_[d] = begin_[d] + n[d];
			b.end_[d] = end_[d] + n[d];
		}
		return b;
	}
	constexpr Box shift(Integer d, Integer n) const {
		std::array<Integer, D> s;
		for (Integer i = 0; i < D; ++i) {
			s[i] = (i == d) ? n : 0_I;
		}
		return shift(s);
	}
	constexpr Box shrink(Integer n) const {
		return expand(-n);
	}
	constexpr Box transpose(Integer d1, Integer d2) const {
		Box t = *this;
		std::swap(t.begin_[d1], t.begin_[d2]);
		std::swap(t.end_[d1], t.end_[d2]);
		return t;
	}
	constexpr Integer flatten(std::array<Integer, D> const &idx) const {
		Integer i = 0_I;
		for (Integer d = 0; d < D; d++) {
			i = span(d) * i + (idx[d] - begin_[d]);
		}
		return i;
	}
	friend constexpr Box bounding(Box const &A, Box const &B) {
		Box I;
		for (Integer d = 0; d < D; ++d) {
			I.begin_[d] = std::min(A.begin_[d], B.begin_[d]);
			I.end_[d] = std::max(A.end_[d], B.end_[d]);
		}
		return I;
	}
	friend constexpr Box intersection(Box const &A, Box const &B) {
		Box I;
		for (Integer d = 0; d < D; ++d) {
			I.begin_[d] = std::max(A.begin_[d], B.begin_[d]);
			I.end_[d] = std::min(A.end_[d], B.end_[d]);
		}
		return I;
	}
	void serialize(auto &arc, unsigned) {
		arc & begin_;
		arc & end_;
	}
};

template <Integer D, typename F>
constexpr void forEach(Box<D> const &box, F const &foo) {
	std::array<Integer, D> idx;
	auto const &beg = box.begin();
	auto const &end = box.end();
	auto const lambda = [&]<Integer I>(auto const &self) {
		if constexpr (I == D) {
			return foo(idx);
		} else {
			for (idx[I] = beg[I]; idx[I] != end[I]; ++idx[I]) {
				self.template operator()<I + 1>(self);
			}
		}
	};
	lambda.template operator()<0>(lambda);
}

template <Integer D, typename F>
constexpr void forEach(Box<D> const &box1, Box<D> const &box2, F const &foo) {
	std::array<Integer, D> idx1, idx2;
	auto const &beg1 = box1.begin();
	auto const &beg2 = box2.begin();
	auto const &end1 = box1.end();
	auto const lambda = [&]<Integer I>(auto const &self) {
		if constexpr (I == D) {
			return foo(idx1, idx2);
		} else {
			for (idx1[I] = beg1[I], idx2[I] = beg2[I]; idx1[I] != end1[I]; ++idx1[I], ++idx2[I]) {
				self.template operator()<I + 1>(self);
			}
		}
	};
	lambda.template operator()<0>(lambda);
}

template <Integer D>
std::ostream &operator<<(std::ostream &os, Box<D> const &box) {
	os << "Box<" << D << ">{ ";
	for (Integer d = 0; d < D; ++d) {
		os << '[' << box.begin()[d] << ',' << box.end()[d] << ')';
		if (d + 1 != D) {
			os << ", ";
		}
	}
	os << " }";
	return os;
}
#endif /* INCLUDE_BOX_HPP_ */
