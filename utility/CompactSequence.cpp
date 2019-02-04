/*
 * CompactSequence.cpp
 *
 *  Created on: 11Dec.,2018
 *      Author: mc7
 */

#include <stdlib.h>
#include <cmath>
#include <string>
#include <sstream>
#include <exception>
#include "CompactSequence.h"

using namespace std;
extern bool _debugging;

namespace flatbush {

const size_t NucsPerSeqWord(4*sizeof(seqwordtype));
		// 4 nucleotides per byte * number of bytes in the sequence word type (e.g., unsigned int)

CompactSequence::CompactSequence(unsigned int len) : sequenceLength(len) {
	/**
	 * Each nucleotide takes up 2 bits, so there are four nucleotides per byte.
	 * Number of bytes required is therefore ceil(length / 4).
	 * We won't assume that each unsigned char is exactly one byte, but we will only use one byte of it anyway.
	 */
	seq.resize(static_cast<size_t>(std::ceil(sequenceLength / static_cast<float>(NucsPerSeqWord))));
}

CompactSequence::CompactSequence(const std::string& str) {
	// initialisation from std::string
	bool _debugging = false;
	sequenceLength = str.size();
	accommodate(sequenceLength);
	encode(str);
}

CompactSequence::CompactSequence(const CompactSequence& cs) {
	// copy constructor
	sequenceLength = cs.getLength();
	seq = cs.getSequence();
}

CompactSequence::~CompactSequence() {
	// TODO Auto-generated destructor stub
}

void CompactSequence::accommodate(unsigned int n) {
	seq.resize(static_cast<size_t>(std::ceil(n / static_cast<float>(NucsPerSeqWord))));
}

void CompactSequence::asDNA(string& str) const {
	decode(str, DNADecoding);
}

void CompactSequence::asRNA(string& str) const {
	decode(str, RNADecoding);
}

void CompactSequence::decode(string& str, const unsigned char decoding[4]) const {
	bool _debugging(true);
	str.clear();
	seqwordtype theWord;
	for (unsigned int i = 0; i < seq.size(); ++i) {
		theWord = seq.at(i);
		for (unsigned int idx = 0; idx < NucsPerSeqWord; ++idx) {
			if (i*NucsPerSeqWord + idx >= sequenceLength) {
				break;
			}
			str.push_back(decoding[theWord & 3]);
			theWord >>= 2;	// now get the next significant two bits
		}
	}
}

void CompactSequence::encode(const string& str) {
	bool _debugging(false);
	unsigned int theWord;
	unsigned int j;
	for (unsigned int i = 0; i < seq.size(); ++i) { // i indexes the unsigned ints stored in seq
		// do each chunk of nucleotide characters as one unsigned int:
		theWord = 0;
		for (unsigned int idx = 0; idx < NucsPerSeqWord; ++idx) {
			j = i * NucsPerSeqWord + idx;
			if (j >= str.size()) {
				break;
			}
			DEBUG(cout << "i = " << i << ", seq-pos = " << j << endl);
			unsigned char ch = toupper(str.at(j));
			DEBUG(cout << "str[" << j << "] = " << ch << endl);
			if (ch != 'A' && ch != 'C' && ch !='G' && ch != 'T' && ch != 'U') {
				stringstream sstr;
				sstr << "Unreadable character \'" << ch << "\' in input string.\n";
				throw new runtime_error(sstr.str());
			}
			theWord += NucleotideEncoding[static_cast<unsigned int>(ch)] << (idx*2);
		}
		seq[i] = theWord;
	}
}

/**
 * These operators return the char at position i in the sequence, NOT the ith unsigned int of the vector.
 */
unsigned char CompactSequence::operator[](unsigned int i) {
	if (i > sequenceLength) {
		throw new range_error("CompactSequence::operator[](unsigned int i): i is out of range.");
	}
	return DNADecoding[(seq.at(i / NucsPerSeqWord) >> 2*(i % NucsPerSeqWord)) & 3];
	// Get the right unsigned int, then shift it right the correct number of bits, then AND with 0x00000011; finally decode assuming it's DNA.
}
const unsigned char CompactSequence::operator[](unsigned int i) const {
	if (i > sequenceLength) {
		throw new range_error("CompactSequence::operator[](unsigned int i): i is out of range.");
	}
	return DNADecoding[(seq.at(i / NucsPerSeqWord) >> 2*(i % NucsPerSeqWord)) & 3];
	// Get the right unsigned int, then shift it right the correct number of bits, then AND with 0x00000011; finally decode assuming it's DNA.
}

CompactSequence& CompactSequence::operator+=(const CompactSequence& cs) {
	bool _debugging(false);
	DEBUG(cout << "this->getLength() = " << getLength() << endl);
	DEBUG(cout << "cs.getLength() = " << cs.getLength() << endl);
	DEBUG(cout << "concatenating two CompactSequences: " << *this << " and " << cs << endl);
	string left;
	asDNA(left);
	string right;
	cs.asDNA(right);
	sequenceLength += cs.getLength();
	accommodate(sequenceLength);
	DEBUG(cout << "right = " << right << endl);
	left += right;
	DEBUG(cout << "left+=right = " << left << endl);
	DEBUG(cout << "sequence length = " << sequenceLength << endl);
	encode(left);
	return *this;
}

CompactSequence& CompactSequence::operator+=(const std::string& right) {
	sequenceLength += right.size();
	accommodate(sequenceLength);
	string left("");
	asDNA(left);	// convert original to DNA string
	left += right;	// append new (hopefully DNA) string
	encode(left);	// encode the whole lot again
	return *this;

}

void CompactSequence::strm(std::ostream& os) const {
	string str;
	asDNA(str);
	os << str;
}

} /* namespace flatbush */
