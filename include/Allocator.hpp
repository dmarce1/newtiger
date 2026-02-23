/*
 * Allocator.hpp
 *
 *  Created on: Feb 22, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_ALLOCATOR_HPP_
#define INCLUDE_ALLOCATOR_HPP_

#include "Integer.hpp"

#include <cstdlib>
#include <functional>

struct Allocator {
	Allocator(Integer size) {
		currentPointer_ = basePointer_ = malloc(size);
	}
	template <typename T>
	T &get() {
		auto *ptr = new (currentPointer_) T();
		delete_.push_back([ptr]() {
			ptr->~T();
		});
		currentPointer_ += sizeof(T);
		return *ptr;
	}
	virtual ~Allocator() {
		while (delete_.size()) {
			delete_.back()();
			delete_.pop_back();
		}
		free(basePointer_);
	}

private:
	void *basePointer_;
	void *currentPointer_;
	std::vector<std::function<void()>> delete_;
};

#endif /* INCLUDE_ALLOCATOR_HPP_ */
