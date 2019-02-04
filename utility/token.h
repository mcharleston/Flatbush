/*
 * token.h
 *
 *  Created on: 2 Jun 2016
 *      Author: mc7
 */

#ifndef TOKEN_H_
#define TOKEN_H_

#include <string>

namespace parsing {

enum tokenType {
	untypedT = 0,
	charT = 1,
	intT = 2,
	doubleT = 3,
	stringT = 4,
	boolT = 5,
	commentT = 6,
	intRangeT = 7
};

class Token {
private:
	int tokNum;
	int lineNumber = -1;
public:
	virtual ~Token();
	Token() : tokNum(0), lineNumber(-1), type(untypedT), b(false), c(' '), d(0.0), i(0), end(0), start(0), s("") {}

	Token(const Token& t) : tokNum(t.tokNum), lineNumber(t.lineNumber), type(t.type), b(t.b), c(t.c), d(t.d),
			i(t.i), end(t.end), start(t.start), s(t.s) {}

	inline int getLineNumber() const { return lineNumber; }

	tokenType getType() const { return type; }

	Token& operator=(const Token& tok);

	friend std::ostream& operator<<(std::ostream & os, const Token &tok);

	inline void set(const std::string & str) { s = str; type = stringT; }
	inline void set(char ch) { c = ch; type = charT; }
	inline void set(int x) { i = x; type = intT; }
	inline void set(int x, int y) { start = x; end = y; type = intRangeT; }
	inline void set(double x) { d = x; type = doubleT; }
	inline void set(bool x) { b = x; type = boolT; }
	inline void setLineNum(int n) { lineNumber = n; }
	inline void setNum(int n) { tokNum = n; }
	std::string toString() const;

	tokenType type;
	bool b;
	char c;
	double d;
	int i;
	int end;
	int start;
	std::string s;

};

} /* namespace flatbush */

#endif /* TOKEN_H_ */
