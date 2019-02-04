/*
 * tokenlist.h
 *
 *  Created on: 2 Jun 2016
 *      Author: mc7
 */

#ifndef TOKENLIST_H_
#define TOKENLIST_H_

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "token.h"

namespace parsing {

//typedef std::vector<parsing::Token>::const_iterator ParseIter;

class TokenList {
private:
	unsigned int lineNumber;
	char ignoreLineChar = '%';
	char beginCommentChar = '[';
	char endCommentChar = ']';
	std::string idchars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_";
public:
	std::vector<parsing::Token> toks;
	std::vector<parsing::Token>::iterator current;

	virtual ~TokenList();
	TokenList() : lineNumber(0) {}
	TokenList(int argn, char** argc);

	void advance() { ++current; }

	inline std::vector<parsing::Token>::const_iterator begin() const { return toks.begin(); }
	inline std::vector<parsing::Token>::iterator begin() { return toks.begin(); }
	inline std::vector<parsing::Token>::const_iterator end() const { return toks.end(); }
	inline std::vector<parsing::Token>::iterator end() { return toks.end(); }

	inline bool isBool() const { return (current->getType() == parsing::boolT); }
	inline bool isChar() const { return (current->getType() == parsing::charT); }
	inline bool isComment() const { return (current->getType() == parsing::commentT); }
	inline bool isInt() const { return (current->getType() == parsing::intT); }
	inline bool isDouble() const { return (current->getType() == parsing::doubleT); }
	inline bool isString() const { return (current->getType() == parsing::stringT); }
	inline bool isIntRange() const { return (current->getType() == parsing::intRangeT); }

	inline char getBeginCommentChar() const { return beginCommentChar; }
	inline char getEndCommentChar() const { return endCommentChar; }

	inline bool hasNext() const { return (current+1 != toks.end()); }
	inline bool hasNextBool() const { return (hasNext() && (current+1)->type == parsing::boolT); }
	inline bool hasNextChar() const { return (hasNext() && (current+1)->type == parsing::charT); }
	inline bool hasNextInt() const { return (hasNext() && (current+1)->type == parsing::intT); }
	inline bool hasNextDouble() const { return (hasNext() && (current+1)->type == parsing::doubleT); }

	inline int numTokens() const { return toks.size(); }

	void reset();

	inline void setBeginCommentChar(char ch) { beginCommentChar = ch; }
	inline void setEndCommentChar(char ch) { endCommentChar = ch; }

	void tokenize(const std::string& fname);
	void tokenize(std::istream& fin);
	void tokenize(std::vector<std::string>& args);
};

} // namespace parsing

//	public void tokenize(String [] args) throws IOException {
//		/**
//		 * This is for tokenizing command-line data; it might be useful elsewhere too
//		 * XXX This is flawed: won't correctly tokenize as args comes from the command-line as
//		 * a very simply constructed array of simple strings separated by spaces.  Duh.
//		 */
//		//		BufferedReader bin = new BufferedReader
//		StringBuilder accum = new StringBuilder();
//		int i;
//		for (i = 0; i < args.length-1; i++) {
//			accum.append(args[i]);
//			accum.append(' ');
//		}
//		accum.append(args[i]);
//		//		accum.append((int) 65535);	// to finish off the reader for the loop in tokenize(BufferedReader)
//		String str = accum.toString();
//		ByteArrayInputStream bins = new ByteArrayInputStream(str.getBytes());
//		InputStreamReader insr = new InputStreamReader(bins);
//		BufferedReader buf = new BufferedReader(insr);
//		tokenize(buf);
//	}



//		/**
//		 * Method to read a number to keep the tokeniser code clearer
//		 *
//		 * @author Lynden Shields
//		 */
//		private char readNumber(StringBuilder w, Token t, Character c, BufferedReader bin) throws IOException {
//			w.append(c);
//			c=readChar(bin);
//			while(Character.isDigit(c)){//Grab the integer part
//				w.append(c);
//				c=readChar(bin);
//			}
//			if(c=='-'){//Check for an integer range
//				bin.mark(1);
//				char d = readChar(bin);
//				bin.reset();
//				if(Character.isDigit(d)){
//					//ignore the hyphen
//					c=readChar(bin);
//					int start = Integer.parseInt(w.toString());
//					w = new StringBuilder();
//					while(Character.isDigit(c)){
//						w.append(c);
//						c= readChar(bin);
//					}
//					t.type = Token.intRangeT;
//					t.set(start,Integer.parseInt(w.toString()));
//					return c;
//				}
//			}
//
//			else if(c=='.'||c=='e'||c=='E'){//if its a double
//				if (c=='.'){//Make sure its not just a full-stop instead of a point
//					bin.mark(1);
//					char d=readChar(bin);
//					bin.reset();
//					if(!Character.isDigit(d)){
//						t.type=Token.intT;
//						t.set(Integer.parseInt(w.toString()));
//						return c;
//					}
//				}
//				w.append(c);
//				c=readChar(bin);
//				while(Character.isDigit(c)){
//					w.append(c);
//					c=readChar(bin);
//				}
//				if(c=='e' || c=='E'){
//					bin.mark(2);			//Read ahead and make sure its an exponent, not the start of a string.
//					char d = readChar(bin);
//					if(Character.isDigit(d)||d=='-'){
//						if(d=='-'){						//Make sure is an int
//							d=readChar(bin);
//							if(!Character.isDigit(d)){ 		//E.G. <INT>"e-"<String>
//								bin.reset();
//								t.type=Token.intT;
//								t.set(Integer.parseInt(w.toString()));
//								return c;
//							}
//						}
//						w.append(c);
//						c= readChar(bin);
//						while(Character.isDigit(c)){
//							w.append(c);
//							c=readChar(bin);
//						}
//					}else{ 									//E.G. <INT>"E"<String>
//						bin.reset();
//						t.type=Token.intT;
//						t.set(Integer.parseInt(w.toString()));
//						return c;
//					}
//					return c;
//				}
//				t.type=Token.doubleT;
//				t.set(Double.parseDouble(w.toString()));
//
//			}else{
//				t.type = Token.intT;
//				t.set(Integer.parseInt(w.toString()));
//			}
//			return c;
//		}
//
//
//		/**
//		 * @param bin
//		 * @return
//		 */
//		private char readChar(BufferedReader bin) {
//			char c = ' ';
//			try {
//				c = (char) bin.read();
//				if (c == '\n' || c == '\r') {
//					lineNum++;
//				}
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//			return c;
//		}
//		public boolean isStringDelimeter(char c) {
//			return (c == '\"' || c == '\'');
//		}
//		public boolean isStringLetter(char c) {
//			if (Character.isLetter(c)) return true;
//			if (c == '_') return true;
//			if (c == '/') return true;
//			return false;
//		}
//		public boolean isWhitespace(char c) {
//			if (c == '\n' || c == '\r' || c == '\t' || c == ' ')
//				return true;
//			return false;
//		}
//		public String toString() {
//			String s = "";
//			for (Token t : toks) {
//				s += t.toString() + "\n";
//			}
//			return s;
//		}
		//	public static void main(String [] args) {
		//		/**
		//		 * This is just to test the token list
		//		 *
		//		 */
		//		TokenList TL = new TokenList();
		//		TL.tokenize("test");
		//
		//	}
//	}
//
//} /* namespace flatbush */

#endif /* TOKENLIST_H_ */
