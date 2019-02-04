/*
 * FASTAParser.cpp
 *
 *  Created on: 20 Jul 2016
 *      Author: mac
 */

#include <algorithm>
#include <initializer_list>
#include <fstream>
#include <sstream>
#include "../utility/debugging.h"
#include "../src/project.h"
#include "FASTAParser.h"

using namespace std;

extern bool _debugging;

namespace parsing {

FASTAParser::FASTAParser(const std::string& fileName, flatbush::Project* pr)
		: proj(pr) {
	ifstream in(fileName, std::ifstream::in);
	/**
	 * tokens are as follows:
	 * '>'
	 * <seqID>
	 * <sequence>
	 *
	 * And grammar is
	 *
	 * <fastafile>     ::= { <fastaSequence> }
	 * <fastaSequence> ::= '>' <seqID> { <seq> }
	 * <seqID>         ::= rest of line up to carriage return
	 * <seq>           ::= entire line unless it begins with '>'
	 */
//	bool _debugging = true;
	string seq, w;
	Token t;
	char c;
	unsigned int lineNum = 1;
	unsigned int tokNum = 0;
	DEBUG(cout << "Tokenizing..." << endl);
	in.get(c);
	while (!in.eof()) {
		// check new lines:
		if (c == '\n' || c == '\r') {
			++lineNum;
		}
		if (c == '>') {
			t.set(c);
			t.type = parsing::charT;
			t.setLineNum(lineNum);
			t.setNum(tokNum);
			TL.toks.push_back(t);
			++tokNum;
			DEBUG(cout << t << endl);
			// skip whitespace:
			in.get(c);
			while (!in.eof() && (c == ' ' || c == '\t')) {
				in.get(c);
			}
			in.unget();
			getline(in, t.s);	// this is the sequence name
			t.setNum(tokNum);
			t.type = parsing::stringT;
			TL.toks.push_back(t);
			++tokNum;
			DEBUG(cout << t << endl);
		} else {
			// This looks innocent enough but there are MANY copy operations!
//			seq = c;	// the first character should be included!
//			getline(in, w);			// here we're putting the line into the string w
//			seq += w;					// now we copy w into seq. That's copy 1.
//			t.set(seq);					// now copy the string into t.s
//			TL.toks.push_back(t);	// now copy it AGAIN. That's THREE copy operations.
//			++tokNum;
//			++lineNum;
			in.unget();
			getline(in, t.s);			// here we're putting the line directly into t.s
			t.setNum(tokNum);
			t.setLineNum(lineNum);
			TL.toks.push_back(t);	// now copy it ONCE. This is acceptable because toks is an STL container and therefore is handling all the resource allocation.
			++tokNum;
			DEBUG(cout << t << endl);
			++lineNum;
		}
		while (!in.eof() && (c == '\n' || c == '\r' || c == '\t' || c == ' ')) {
			if (c == '\n' || c == '\r') {
				++lineNum;
			}
			in.get(c);
			DEBUG(cout << "reading char " << c << endl);
		}
		w = "";
		in.get(c);
	}
}

void FASTAParser::parse() {
	/*
	 * FASTA grammar:
	 * <fasta> ::= { <item> }
	 * <item> :: '>' <id> <sequence>
	 */
	bool _debugging = false;
	DEBUG(cout << "FASTAParser::parse()" << endl);
	TL.reset();
	while (hasNext()) {
		eat('>');				// >
		string id = getString();	// sequenceName, should be complete line
		DEBUG(cout << "seqID = " << id << endl);
		advance();
		stringstream seq;
		seq << getString();	// there must be at least one line!
		while (!matches('>')) {
			if (hasNext()) {
				advance();
			} else {
				break;
			}
			if (matches('>')) {
				break;
			}
			seq << getString();
			DEBUG(cout << "complete sequence is " << seq.str() << endl);
		}
//		alignment.insert(pair<string, string>(id, seq.str()));
		string str(seq.str());
		std::transform(str.begin(), str.end(), str.begin(), ::toupper);	// have to convert to uppercase. sheesh.
		proj->addRawSequence(id, str);
	}
}

//std::initializer_list<std::string> FASTAParser::getSuffixes() {
//	return { ".fa", ".fst", ".fasta" };
//}

} /* namespace parsing */
