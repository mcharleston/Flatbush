/*
 * Node.cpp
 *
 *  Created on: 18 Apr 2016
 *      Author: mc7
 */

#include "node.h"
#include "parser.h"

using namespace std;

namespace flatbush {

std::string Node::labelChars("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_");
std::string Node::alpha("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");

Node::~Node() { }

Node::Node() : label(), parent(nullptr), firstChild(nullptr), sibling(nullptr),
		branchLength(0.0), height(1), idx(0), _marked(false) { }

Node::Node(const string& str) : label(str), parent(nullptr), firstChild(nullptr), sibling(nullptr),
		branchLength(0.0), height(1), idx(0), _marked(false) { }

void Node::addChild(Node* sib) {
	sib->setParent(this);
	if (firstChild == nullptr) {
		firstChild = sib;
	} else {
		Node* child = firstChild;
		while (child->sibling != nullptr) {
			child = child->sibling;
		}
		child->sibling = sib;
	}
}

int Node::calculateHeight() {
	height = 0;
	if (isLeaf()) {
		height = 0;
	} else {
		Node* child = firstChild;
		while (child != nullptr) {
			height = max(height, child->getHeight()+1);
			child = child->getSibling();
		}
	}
	return height;
}

void Node::clearTraversal() {
	// unmark all nodes.
	_marked = false;
	Node* child = firstChild;
	while (child != nullptr) {
		child->clearTraversal();
		child = child->sibling;
	}
}

set<Node*>* Node::getAllDescendants(set<leafset>& splits, unsigned int numLeaves) {
	if (isLeaf()) {
		if (!_marked) {
			leafset l(numLeaves);
			descendantLeaves.insert(this);
			l[this->idx] = true;
			cout << "{ " << this->idx << " }" << endl;
			splits.insert(l);
			_marked = true;
		}
		return &descendantLeaves;
	}
	descendantLeaves.clear();
	Node* child = firstChild;
	leafset l(numLeaves);
	while (child != nullptr) {
		set<Node*>* desc = child->getAllDescendants(splits, numLeaves);
		descendantLeaves.insert(desc->begin(), desc->end());
		for (Node* n : child->descendantLeaves) {
			l[n->idx] = true;
		}
		splits.insert(l);
		child = child->sibling;
	}
	if (!_marked) {
		cout << "{ ";
		for (unsigned int i = 0; i < numLeaves; ++i) {
			if (l[i]) {
				cout << i << " ";
			}
		}
		cout << "}" << endl;
		_marked = true;
	}
	return &descendantLeaves;
}

void Node::setChildren(Node* left, Node* right) {
	firstChild = left;
	firstChild->sibling = right;
	left->parent = this;
	right->parent = this;
}

} // namespace parsing
