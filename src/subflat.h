/*
 * subflat.h
 *
 *  Created on: 22 Jul 2016
 *      Author: mac
 */

#ifndef SRC_SUBFLAT_H_
#define SRC_SUBFLAT_H_

#include <map>
#include <vector>
#include "../utility/pattern.h"
//#include "../utility/alignment.h"
//#include "../utility/node.h"

namespace flatbush {

class Node;

typedef std::vector<bool> leafset;

unsigned int Fitch(const unsigned int* T, const flatbush::leafset & C, unsigned int n);

std::string splash();
//void putHalfSplitMask(std::vector<Pattern<char>>& idx, flatbush::leafset X);
std::vector<double> svdSplit(std::map<Pattern<unsigned int>, int> signedCount, unsigned int n, std::vector<bool> A);
void putAllNontrivialHalfSplits(std::vector<flatbush::leafset>& S, int ntax, bool _verbose);
void readNewickTree(flatbush::Node * n, const char* treestr);
void timeSVD();

}

#endif /* SRC_SUBFLAT_H_ */

