/*
 * parser.cpp
 *
 *  Created on: 2 Jun 2016
 *      Author: mc7
 */


#include "parser.h"
#include "../src/project.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <set>
#include <sstream>

#include "debugging.h"
#include "appexception.h"
#include "token.h"

using namespace std;
using namespace parsing;

extern bool _debugging;

Parser::Parser() 	: _ignoreCase(true), _verbose(false), tokNumCount(0), proj(nullptr) {
	DEBUG(cout << "instantiating Parser object" << endl);
}

Parser::Parser(const string& filename) : _ignoreCase(true), _verbose(false), tokNumCount(0), proj(nullptr) {
	ifstream infile(filename);
	if (!infile.is_open()) {
		throw new flatbush::app_exception("Cannot open this file: are you sure it exists?");
	}
	TL.tokenize(infile);
}

// Move along the token list to the next one, if there is one.
void Parser::advance() {
	bool _debugging = false;
	if (current().getType() == parsing::tokenType::commentT) {
		while (TL.hasNext() && current().getType() == parsing::tokenType::commentT) {
			TL.advance();
		}
	} else {
		if (getTL().hasNext()) {
			getTL().advance();
		}
	}
	DEBUG (cout << current() << endl);
}

bool Parser::assignNextBoolean() {
	advance();
	ignore('=');
	return getBoolean();
}
double Parser::assignNextDouble() {
	advance();
	ignore('=');
	return getDouble();
}
int Parser::assignNextInt() {
	advance();
	ignore('=');
	return getInt();
}
std::string Parser::assignNextString() {
	advance();
	ignore('=');
	return getString();
}
std::vector<double> Parser::assignNextVector() {
	advance();
	ignore('=');
	return getVector();
}

bool Parser::currentIs(char ch) {
	if (current().type == parsing::tokenType::charT && current().c == ch)
		return true;
	string w;
	w += ch;
	if (current().type == parsing::tokenType::stringT &&
			eq(w, current().s))
		return true;
	return false;
}
bool Parser::currentIs(int i) {
	if (current().type == parsing::tokenType::intT && current().i == i) {
			return true;
	}
	return false;
}
bool Parser::currentIs(std::string s) {
	if (current().type == parsing::tokenType::stringT && eq(current().s, s))
		return true;
	if (s.length() == 1 && current().type == parsing::tokenType::charT) {
		if (_ignoreCase) {
			return (::tolower(current().s[0]) == ::tolower(s[0]));
		} else {
			return (current().s[0] == s[0]);
		}
	}
	return false;
}
bool Parser::currentIs(std::initializer_list<std::string> ss) {
	if (!TL.hasNext())
		return false;
	for (std::string s : ss) {
		if (currentIs(s))
			return true;
	}
	return false;
}

void Parser::eat(char ch) {
	expect(ch);
	advance();
}
void Parser::eat(std::string s) {
	expect(s); // this is a required string
	advance();
}
void Parser::eat(std::initializer_list<std::string> ss) {
	for (std::string s : ss) {
		if (currentIs(s)) {
			advance();
			return;
		}
	}
	std::string msg("Parsing error: \"");
	msg += current().s;
	msg += "\" not found in list of acceptable options";
	throw new parse_error(msg);
}

bool Parser::eq(string str1, string str2) {
	if (_ignoreCase) {
		transform(str1.begin(), str1.end(), str1.begin(), ::tolower);
		transform(str2.begin(), str2.end(), str2.begin(), ::tolower);
	}
	return (str1.compare(str2) == 0);
}

void Parser::expect(char ch) {
	skipComments();
	if (!currentIs(ch)) {
		stringstream msg;
		msg << "On line #" << TL.current->getLineNumber() << ", expecting character \'"
				<< ch << "\' here but got token #" << TL.current->toString() << ".";
		throw new parse_error(msg.str());
	}
}
void Parser::expect(std::string s) {
	skipComments();
	if (!eq(s, TL.current->s)) {
		stringstream msg;
		msg << "On line #" << TL.current->getLineNumber() << ", expecting string \""
				<< s << "\" here but got token #" << TL.current->toString() << ".";
		throw new parse_error(msg.str());
	}
}

bool Parser::getBoolean() {
	skipComments();
	ignore('=');
	if (currentIs(1) || currentIs('t') || currentIs('y') || matches({"true", "yes"})) {
		return true;
	}
	return false;
}
char Parser::getChar() {
	skipComments();
	if (current().type == parsing::tokenType::charT) {
		return current().c;
	}
	if (parsing::tokenType::stringT == current().type && current().s.length() == 1) {
		return current().s[0];
	}
	std::string msg("Expecting a char here but got nothing I could use: ");
	msg += TL.current->toString();
	throw new parse_error(msg);
}
double Parser::getDouble() {
	double d;
	bool _debugging(true);
	skipComments();
	DEBUG(cout << "the token is " << current().toString() << endl);
	if (current().type == parsing::tokenType::doubleT) {
		d = current().d;
	} else if (current().type == parsing::tokenType::intT) {
		d = (double) current().i;
	} else {
		std::stringstream ss;
		ss << "On input line " << current().getLineNumber() << ", expecting a double here but got nothing I could use: "
				<< TL.current->toString();
		throw new parse_error(ss.str());
	}
	return d;
}
int Parser::getInt() {
	skipComments();
	if (current().type != parsing::tokenType::intT) {
		std::string msg("Expecting an int here but got nothing I could use: ");
		msg += TL.current->toString();
		throw new parse_error(msg);
	}
	int i = current().i;
	return i;
}
std::string Parser::getString() {
	std::string s = "";
	skipComments();
	if (current().type == parsing::tokenType::stringT) {
		s = current().s;
		return s;
	} else if (current().type == parsing::tokenType::charT) {
		s += current().c;
		return s;
	} else if (current().type == parsing::tokenType::intT) {
		s += std::to_string(current().i);
		return s;
	} else if (current().type == parsing::tokenType::doubleT) {
		s += std::to_string(current().d);
		return s;
	} else {
		std::string msg("Expecting a string here but got nothing I could use: ");
		msg += TL.current->toString();
		throw new parse_error(msg);
	}
	return s;
}
std::vector<double> Parser::getVector() {
	skipComments();
	eat('(');
	std::vector<double> vec;
	while (!currentIs(')')) {
		vec.push_back(getDouble());
		advance();
	}
	eat(')');
	return vec;
}

bool Parser::hasNext(char ch) {
	if (!hasNext())
		return false;
	if (current().type == parsing::tokenType::charT && current().c == ch)
		return true;
	return false;
}
bool Parser::hasNext(double d) {
	if (!hasNext())
		return false;
	if (current().type == parsing::tokenType::doubleT && current().d == d)
		return true;
	return false;
}
bool Parser::hasNext(int i) {
	if (!hasNext())
		return false;
	if (current().type == parsing::tokenType::intT && current().i == i)
		return true;
	return false;
}

void Parser::ignore(char ch) {
	if (hasNext(ch))
		advance();
}
void Parser::ignore(std::string s) {
	if (currentIs(s))
		advance();
}

void Parser::ignore(int n_args, std::initializer_list<std::string> ss) {
	for (std::string s : ss) {
		if (currentIs(s)) {
			advance();
			return;
		}
	}
}

bool Parser::matches(const char c) {
	return currentIs(c);
}
bool Parser::matches(const string& s) {
	return currentIs(s);
}
bool Parser::matches(const initializer_list<char>& cc) {
	for (char c : cc) {
		if (currentIs(c)) {
			return true;
		}
	}
	return false;
}
bool Parser::matches(const initializer_list<string>& ss) {
	for (std::string s : ss) {
		if (currentIs(s)) {
			return true;
		}
	}
	return false;
}

//double Parser::nextBool() {
//	return 0.0;
//}
//
//std::string Parser::nextString() {
//	advance();
//	return getString();
//}
//
//const string& Parser::nextToken(vector<string> &tokens, parseIter iter) {
//	DEBUG(cout << "token:" << **iter << "; next is ");
//	++(*iter);
//	DEBUG(cout << **iter << endl);
//	skipComments(tokens, iter);
//	return **iter;
//}
//
//void Parser::skipTo(std::string s) {
//	while (current().type != parsing::stringT ||
//			(current().type == parsing::stringT && !eq(current().s, s))) {
//		if (hasNext())
//			advance();
//		else
//			break;
//	}
//}

void Parser::skipComments() {
	bool _debugging = false;
	DEBUG(cout << "Skipping comments...\n");
	int n = 0;
	DEBUG(cout << "this token: " << current() << endl);
	while (current().getType() == commentT) {
		DEBUG(cout << "Comment identified: " << current() << endl);
		++n;
		advance();
		DEBUG(cout << "next token: " << current() << endl);
	}
	DEBUG(cout << n << " comments skipped.\n");
}

void Parser::skipWhitespace() {
	// hah. Nothing here.
}

bool matches(char ch, std::initializer_list<char>& cc) {
	for (char c : cc) {
		if (ch == c) {
			return true;
		}
	}
	return false;
}
bool matches(const string&str, const std::initializer_list<std::string>& patterns) {
	vector<string> vec(patterns.begin(), patterns.end());
//	return matches(str, vec);
	for (std::string s : patterns) {
		if (s.compare(str) == 0) {
			return true;
		}
	}
	return false;
}
//bool matches(const string&str, const std::vector<std::string>& patterns) {
//}

bool matchesIgnoreCase(char ch, std::initializer_list<char>& cc) {
	ch = tolower(ch);
	for (char c : cc) {
		if (ch == tolower(c)) {
			return true;
		}
	}
	return false;
}
bool matchesIgnoreCase(string str, initializer_list<string> patterns) {
	// Convert initializer_list to vector and call the (string, vector) overload.
	transform(str.begin(), str.end(), str.begin(), ::tolower);
	vector<string> vec(patterns.begin(), patterns.end());
//	return matchesIgnoreCase(str, vec);
	for (string s : patterns) {
		transform(s.begin(), s.end(), s.begin(), ::tolower);
		if (str.compare(s) == 0) {
			return true;
		}
	}
	return false;
}

//ostream& operator<<(ostream& os, const set<string>& S) {
//	os << "{ ";
//	for (string s : S) {
//		os << s << " ";
//	}
//	os << "}";
//	return os;
//}
