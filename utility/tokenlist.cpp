/*
 * tokenlist.cpp
 *
 *  Created on: 2 Jun 2016
 *      Author: mc7
 */

#include "debugging.h"
#include "tokenlist.h"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace std;

namespace parsing {

TokenList::~TokenList() {
	// TODO Auto-generated destructor stub
}

TokenList::TokenList(int argn, char** argc) {
	stringstream ss;
	for (int i = 1; i < argn; ++i) {
		ss << argc[i] << " ";
	}
	tokenize(ss);
}

void TokenList::reset() {
	current = toks.begin();
}


void TokenList::tokenize(const string& fname) {
	/**
	 * Read the text in a file and turn it into a list of tokens
	 *
	 * Read each character in turn and then allocate to
	 *   number
	 *   string
	 *   char (of different types)
	 */
	// NOT AT ALL EXCEPTION SAFE
	ifstream fin(fname);
	tokenize(fin);
}

void TokenList::tokenize(istream& in) {
	bool _debugging = true;
	string w;
	Token t;
	char c;
	unsigned int lineNum = 1;
	unsigned int tokNum = 0;
	in.get(c);
	while (!in.eof()) {
		// check new lines:
		if (c == '\n' || c == '\r') {
			++lineNum;
		}
		// check for character that says "ignore this line"
		if (c == ignoreLineChar) {
			t.type = parsing::commentT;
			getline(in, w);
			t.set(w);
			++lineNum;
		} else if (c == beginCommentChar) {
			// begin comment: keep going until we get to the end comment char:
			w += c;
			unsigned int commentDepth = 1;
			// allow nested comments:
			while (!in.eof() && commentDepth > 0) {
				in.get(c);
				w += c;
				if (c == beginCommentChar) {
					++commentDepth;
				} else if (c == endCommentChar) {
					--commentDepth;
				}
			}
			t.set(w);
			t.type = parsing::commentT;
			in.get(c);
		} else if (isdigit(c) || c == '-') {
			// try to read this in as a number
			bool _containsDecimalPoint(false);
			while (!in.eof() && (isdigit(c) || c == '.' || c == '-')) {
				w += c;
				if (c == '.') {
					_containsDecimalPoint = true;
				}
				in.get(c);
			}
			try {
				if (_containsDecimalPoint) {
					t.set(std::stod(w));
				} else {
					t.set(std::stoi(w));
				}
			} catch (const std::invalid_argument& e) {
				// we'll try for an int next
				try { // float?
					t.set(std::stod(w));
				} catch (const std::invalid_argument& e) {
					// ok, treat it as a string
					t.set(w);
				}
			}
//			in.get(c);
		} else if (c == '\"') {
			// quote-delimited string
			in.get(c);
			while (!in.eof() && (c != '\"')) {
				w += c;
				in.get(c);
			}
			t.type = parsing::stringT;
			t.set(w);
			in.get(c);
		} else if (idchars.find(c) != string::npos) {
			// ordinary string
			while (!in.eof() && idchars.find(c) != string::npos) { //(c != ' ' && c != '\t')
				w += c;
				in.get(c);
			}
			t.type = parsing::stringT;
			t.set(w);
		} else {
			// special character
			w = c;
			t.type = parsing::charT;
			t.set(c);
			in.get(c);
		}
		t.setLineNum(lineNum);
		t.setNum(tokNum);
		toks.push_back(t);
		while (!in.eof() && (c == '\n' || c == '\r' || c == '\t' || c == ' ')) {
			if (c == '\n' || c == '\r') {
				++lineNumber;
			}
			in.get(c);
		}
		w = "";
		++tokNum;
		DEBUG(cout << t << endl);
	}
}

} /* namespace flatbush */
