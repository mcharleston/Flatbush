/*
 * parser.h
 *
 *  Created on: 2 Jun 2016
 *      Author: mc7
 */

#ifndef PARSER_H_
#define PARSER_H_

#include <initializer_list>
#include <iostream>
#include <iterator>
#include <set>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <vector>

#include "tokenlist.h"

namespace flatbush {
class Project;
}

namespace parsing {

class parse_error : public std::runtime_error {
public:
	explicit parse_error(const std::string& what_arg) : std::runtime_error(what_arg) {
		std::cerr << what_arg << std::endl;
	}
	explicit parse_error(const char* what_arg) : std::runtime_error(what_arg) { }
};

class Parser {
protected:
	bool _ignoreCase;
	bool _verbose;
	TokenList TL;
	unsigned int tokNumCount;
	flatbush::Project* proj;
public:
	virtual ~Parser() {}
	explicit Parser();
	explicit Parser(const std::string& fileName);
	explicit Parser(TokenList& tl) : _ignoreCase(true), _verbose(false), tokNumCount(0), proj(nullptr) { TL = tl; }

	void advance();

	bool assignNextBoolean();
	double assignNextDouble();
	int assignNextInt();
	std::string assignNextString();
	std::vector<double> assignNextVector();

	inline Token& current() { return *(TL.current); }

	bool currentIs(char ch);
	bool currentIs(int i);
	bool currentIs(std::string s);
	bool currentIs(std::vector<std::string> ss);
	bool currentIs(std::initializer_list<std::string> ss);

	void eat(char ch);
	void eat(std::string s);
	void eat(std::initializer_list<std::string> ss);

	bool eq(std::string str1, std::string str2);

	void expect(char ch);
	void expect(std::string s);

	bool getBoolean();
	char getChar();
	double getDouble();
	int getInt();
//	std::string getLine() {
//		// read a set of characters to the end of the line. ignore gaps.
//		std::string s;
//		while (!currentIs("\n")) {
//			s += nextString();
//		}
//		return s;
//	}

	std::string getString();
	inline TokenList& getTL() { return TL; }
	std::vector<double> getVector();

	inline bool hasNext() { return TL.hasNext(); }
	bool hasNext(char ch);
	bool hasNext(double d);
	bool hasNext(int i);
//	std::string nextString();

	void ignore(char ch);
	void ignore(std::string s);
	void ignore(int n_args, std::initializer_list<std::string> ss);

	inline bool isBool() const { return TL.isBool(); }
	inline bool isChar() const { return TL.isChar(); }
	inline bool isComment() const { return TL.isComment(); }
	inline bool isDouble() const { return TL.isDouble(); }
	inline bool isInt() const { return TL.isInt(); }
	inline bool isNumber() { return (isDouble() || isInt() ); }
	inline bool isString() const { return TL.isString(); }

	bool matches(const char c);
	bool matches(const std::string& s);
	bool matches(const std::initializer_list<char>& cc);
	bool matches(const std::initializer_list<std::string>& ss);

	inline Token& next() { return *(TL.current+1); }

	virtual void parse() {};
	void setProject(flatbush::Project* p) { proj = p; }

	void setVerbose(bool b) { _verbose = b; }
	void skipComments();
	void skipTo(std::string s);
	void skipWhitespace();

	};

void toLower(std::string& str);

bool matches(char ch, std::initializer_list<char>& chars);
bool matches(const std::string& str, const std::initializer_list<std::string>& patterns);
//bool matches(const std::string& str, const std::vector<std::string>& patterns);

bool matchesIgnoreCase(char ch, std::initializer_list<char>& chars);
bool matchesIgnoreCase(std::string str, std::initializer_list<std::string> patterns);
//bool matchesIgnoreCase(std::string str, std::vector<std::string>& patterns);

//std::ostream& operator<<(std::ostream& os, const std::set<std::string>& S);
}
#endif /* PARSER_H_ */
