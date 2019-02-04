/*
 * FASTAParser.h
 *
 *  Created on: 20 Jul 2016
 *      Author: mac
 */

#ifndef UTILITY_FASTAPARSER_H_
#define UTILITY_FASTAPARSER_H_

#include <initializer_list>
#include <map>
#include <vector>
#include "../src/project.h"
#include "parser.h"

namespace parsing {

const std::vector<std::string> FASTASuffixes = { ".fa", ".fas", ".fst", ".fasta" };

class FASTAParser : public Parser {
private:
//	std::map<std::string, std::string> alignment;
	flatbush::Project* proj;
public:
	virtual ~FASTAParser() {}
	FASTAParser() : proj(nullptr) {}
	FASTAParser(const std::string& fileName, flatbush::Project* p);
	FASTAParser(TokenList& tl) : Parser(tl), proj(nullptr) {}

//	inline std::map<std::string, std::string> getAlignment() { return alignment; }

	void parse();
	// put the acceptable FASTA format suffixes in here as a vector
//	static std::initializer_list<std::string> getSuffixes();
};

} /* namespace parsing */

#endif /* UTILITY_FASTAPARSER_H_ */
