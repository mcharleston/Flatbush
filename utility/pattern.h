/*
 * pattern.h
 *
 *  Created on: 24 May 2016
 *      Author: mc7
 */

#ifndef PATTERN_H_
#define PATTERN_H_

#include <iostream>
#include <sstream>
#include <vector>

template <typename T>
class Pattern {
private:
	static const unsigned int seed = 1;
	std::vector<T>* tuple;
public:
	// destructor:
	~Pattern() {
		delete tuple;
	}

	// default constructor:
	Pattern() {
		tuple = new std::vector<T>();
	}

	// normal constructor:
	Pattern(unsigned int n, T t = 0) {
		tuple = new std::vector<T>(n, t);
	}

	// normal constructor:
	Pattern(const std::string& pat) {
		tuple = new std::vector<T>(pat.begin(), pat.end());
	}

	// copy constructor:
	Pattern(const Pattern<T>& pat) {
		tuple = new std::vector<T>(*(pat.tuple));
	}

	// constructor from initializer list:
	Pattern(std::initializer_list<T> list) {
		tuple = new std::vector<T>();
		for (T t : list) {
			tuple->push_back(t);
		}
		tuple->shrink_to_fit();
	}

	/**
	 * Functions:
	 */
	size_t size() const { return tuple->size(); }
	inline void append(const T & item) { tuple->push_back(item); }
	inline void clear() { tuple->clear(); }

	/**
	 * Operators:
	 */
	// copy operator
	Pattern<T>& operator=(const Pattern<T> p) {
		tuple = new std::vector<T>(p.tuple);
		return *this;
	}
	Pattern<T>& operator=(T t) {
		for (unsigned int i = 0; i < tuple->size(); ++i) {
			(*tuple)[i] = t;
		}
		return *this;
	}
	// access:
	T& operator[](unsigned int i) { return tuple->at(i); }
	const T& operator[](unsigned int i) const { return tuple->at(i); }
	T& get(unsigned int i) { return tuple->at(i); }

	// comparison:
	bool operator<(const Pattern& pat) const {
		// shorter patterns are "before" longer ones; aside from that, use lex ordering on T
		if (tuple->size() != pat.tuple->size()) {
			return tuple->size() < pat.tuple->size() ? true : false;
		}
		size_t comparisonLength = pat.tuple->size();
		if (pat.tuple->size() < comparisonLength) {
			comparisonLength = pat.tuple->size();
		}
		for (unsigned int i = 0; i < comparisonLength; ++i) {
			if ((*tuple)[i] != (*pat.tuple)[i]) {
				return (*tuple)[i] < (*pat.tuple)[i] ? true : false;
			}
		}
		return false;
	}

	friend std::ostream& operator<<(std::ostream& os, const Pattern<T>& pat) {
		os << "{";
		for (T t : *(pat.tuple)) {
			os << t;
		}
		os <<"}";
		return os;
	}

	void resize(size_t n) {
		tuple->resize(n);
	}

	void strm(std::ostream & os = std::cout) { os << *this; }

};

// hash_code
namespace std {
	template <typename T> struct hash< Pattern<T> > {
		size_t operator()(const Pattern<T> & pat) const {
			std::hash<std::string> str_hash;
			std::stringstream ss;
			ss << pat;
			return str_hash(ss.str());
//			std::hash<T> t_hash;
//			unsigned int h = 1;
//			unsigned int prime = 104729;	// a smallish (the 10,000th) prime, should do for now...
//			size_t len = pat.size();
//			for (size_t i = 0; i < len; ++i) {
//				h ^= prime * h + t_hash(pat[i]);
//			}
//			return h % prime;
		}
	};
}

#endif /* PATTERN_H_ */
