/*
 * token.cpp
 *
 *  Created on: 2 Jun 2016
 *      Author: mc7
 */

#include "../utility/token.h"

#include <iostream>
#include <string>
#include <sstream>

using namespace std;

namespace parsing {

Token::~Token() {}

Token& Token::operator=(const Token& tok) {
	tokNum = tok.tokNum;
	lineNumber = tok.lineNumber;
	type = tok.type;
	b = tok.b;
	c = tok.c;
	d = tok.d;
	i = tok.i;
	start = tok.start;
	end = tok.end;
	s = tok.s;
	return *this;
}

std::string Token::toString() const {
	stringstream str("#");
	str << tokNum << " type=";
	switch (type) {
		case untypedT: str << "untyped"; break;
		case charT: str << "c \'" << c << "\'"; break;
		case intT: str << "i " << i; break;
		case doubleT: str << "d " << d; break;
		case stringT: str << "s \"" << s << "\""; break;
		case boolT: str << "b " << b; break;
		case commentT: str << "comment \"" << s << "\""; break;
		case intRangeT: str << "range " << start << " - " << end; break;
	}
	return str.str();
}


ostream& operator<<(ostream & os, const Token &tok) {
	os << tok.toString();
	return os;
}

} /* namespace flatbush */
