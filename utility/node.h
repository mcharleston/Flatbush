/*
 * Node.h
 *
 *  Created on: 18 Apr 2016
 *      Author: mc7
 */

#ifndef NODE_H_
#define NODE_H_

#include <set>
#include <string>
#include <vector>
#include "tokenlist.h"
#include "../src/subflat.h"

namespace flatbush {

class Node {
private:
	std::string label;
	Node* parent;
	Node* firstChild;
	Node* sibling;
	double branchLength;
	int height;
	int idx;
	std::set<Node*> descendantLeaves;
	bool _marked;
public:
	virtual ~Node();
	Node();
	Node(const std::string& str);

	void addChild(Node* child);

	int calculateHeight();

	void clearTraversal();
	std::set<Node*>* getAllDescendants(std::set<leafset>& splits, unsigned int numLeaves);
	inline Node* getFirstChild() { return firstChild; }
	inline std::string getLabel() { return label; }
	inline int getHeight() { return height; }
	inline unsigned int getIndex() const { return idx; }
	Node* getSibling() { return sibling; }
	inline Node* getParent() { return parent; }

//	void parse(const TokenList& TL);
	inline bool isLeaf() { return (firstChild == nullptr); }

	inline void setBranchLength(double d) { branchLength = d; }
	void setChildren(Node* left, Node* right);
	inline void setIndex(unsigned int i) { idx = i; }
	void setLabel(std::string str) { label = str; }
	inline void setParent(Node* p) { parent = p; }

	static std::string labelChars;
	static std::string alpha;

};

} /* namespace flatbush */

#endif /* NODE_H_ */
