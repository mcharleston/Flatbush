/*
 * NEXUSParser.h
 *
 *  Created on: 20 Jul 2016
 *      Author: mac
 */

#ifndef NEXUSPARSER_H_
#define NEXUSPARSER_H_

//#include "alignment.h"
#include "../src/subflat.h"
#include "node.h"
#include "parser.h"
#include "phylogeny.h"

namespace parsing {

const std::vector<std::string> NEXUSSuffixes = { ".nx", ".nxs", ".nex", ".nexus" };

class NEXUSParser : public Parser {
//	flatbush::Alignment* A;
	flatbush::Project* proj;
public:
	virtual ~NEXUSParser() { }
	NEXUSParser(const std::string& fileName, flatbush::Project* pr);
	NEXUSParser(TokenList& tl, flatbush::Project* pr) : Parser(tl), proj(pr) {}

	void parse();
	void parseBranchLength(flatbush::Node* v);
	void parseDataBlock();
	void parseDataFormat();
	void parseDimensionsComponent();
	void parseFlatbushBlock();
	void parseMatrix();
	void parseNEXUSBlock();
	void parseNewickFormatTree(flatbush::Phylogeny* T);
	void parseNewickSubtree(flatbush::Node* v);
	void parseSetsBlock();
	void parseTaxaBlock();

	void skipBlock();
};

} /* namespace parsing */

#endif /* NEXUSPARSER_H_ */
