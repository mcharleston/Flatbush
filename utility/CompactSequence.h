/*
 * CompactSequence.h
 *
 *  Created on: 11Dec.,2018
 *      Author: mc7
 */

#ifndef UTILITY_COMPACTSEQUENCE_H_
#define UTILITY_COMPACTSEQUENCE_H_

#include <vector>
#include <stdlib.h>
#include <iostream>
#include "debugging.h"

namespace flatbush {

typedef unsigned int seqwordtype;

class CompactSequence {
private:
	size_t sequenceLength;
	std::vector<seqwordtype> seq;
public:
	explicit CompactSequence() : sequenceLength(0) {}			// default constructor is marked private and explicit to avoid accidental creation.
	CompactSequence(unsigned int len);	// empty constructor
	CompactSequence(const std::string& str);	// initialisation from std::string
	CompactSequence(const CompactSequence& cs);	// copy constructor
	virtual ~CompactSequence();
	void accommodate(unsigned int n);
	void asDNA(std::string& str) const;
	void asRNA(std::string& str) const;
	void decode(std::string& str, const unsigned char decoding[4]) const;
	void encode(const std::string& str);
	CompactSequence& operator+=(const CompactSequence& cs);
	CompactSequence& operator+=(const std::string& str);
	inline size_t getLength() const { return sequenceLength; }
	inline std::vector<seqwordtype> const & getSequence() const { return seq; }
	// access:
	unsigned char operator[](unsigned int i);
	const unsigned char operator[](unsigned int i) const;
	seqwordtype at(unsigned int i);	// does range-checking

	// equality:
	bool operator==(const CompactSequence& cs) const { return seq == cs.seq; }
	// inequality:
	bool operator<(const CompactSequence& cs) const {
		return seq < cs.seq;
//		// shorter patterns are "before" longer ones; aside from that, use lex ordering on T
//		if (length != cs.length) {
//			return length < cs.length ? true : false;
//		}
//		size_t comparisonLength = std::min(length, cs.length);
//		for (unsigned int i = 0; i < comparisonLength; ++i) {
//			if ((*tuple)[i] != (*pat.tuple)[i]) {
//				return (*tuple)[i] < (*pat.tuple)[i] ? true : false;
//			}
//		}
//		return false;
	}
	bool operator>(const CompactSequence& cs) const {
		return cs < *this;
	}

//	std::size_t operator()() const {
//		return std::hash<vector<unsigned char>>(seq);
//	}

	friend std::ostream& operator<<(std::ostream& os, const CompactSequence& cs) {
		std::string str;
		cs.asDNA(str);
		os << "{" << str << "}";
		return os;
	}

	void strm(std::ostream& os = std::cout) const;
};

/**
 * ASCII encoding of ACGTU:
 * Character	ASCII value		encoding for compact sequence
 * 	A				65								0
 * 	a				97								0
 * 	C				67								2
 * 	c				99								2
 * 	G				71								1
 * 	g				103							1
 * 	T				84								3
 * 	t				116							3
 * 	U				85								3
 * 	u				117							3
 * 	Everything else:							33
 */

const unsigned char NucleotideEncoding[256] = {
		33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, //   0 -  15
		33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, //  16 -  31
		33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, //  32 -  47
		33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, //  48 -  63
		//   A       C               G
		33,  0, 33,  2, 33, 33, 33,  1, 33, 33, 33, 33, 33, 33, 33, 33, //  64 -  79
		//               T   U
		33, 33, 33, 33,  3,  3, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, //  80 -  95
		//   a       c               g
		33,  0, 33,  2, 33, 33, 33,  1, 33, 33, 33, 33, 33, 33, 33, 33, //  96 - 111
		//               t   u
		33, 33, 33, 33,  3,  3, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, // 112 - 127
		33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, // 128 - 143
		33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, // 144 - 159
		33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, // 160 - 175
		33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, // 176 - 191
		33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, // 192 - 207
		33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, // 208 - 223
		33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, // 224 - 239
		33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33  // 240 - 255
};

const unsigned char DNADecoding[4] = { 'A', 'G', 'C', 'T' };
const unsigned char RNADecoding[4] = { 'A', 'G', 'C', 'U' };

} /* namespace flatbush */

namespace std {
	template<> struct
	hash<flatbush::CompactSequence> {
		std::size_t operator()(const flatbush::CompactSequence & cs) const {
			const std::vector<flatbush::seqwordtype>& seq = cs.getSequence();
			std::size_t seed = seq.size();
			for (auto &i : seq) {
				seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
		}
	};
}


#endif /* UTILITY_COMPACTSEQUENCE_H_ */
