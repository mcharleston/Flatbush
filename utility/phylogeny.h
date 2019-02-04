/*
 * phylogeny.h
 *
 *  Created on: 20 Jul 2016
 *      Author: mac
 */

#ifndef UTILITY_PHYLOGENY_H_
#define UTILITY_PHYLOGENY_H_

#include <map>
#include <stdio.h>
#include <vector>
#include "../src/project.h"
#include "node.h"

namespace flatbush {

class Phylogeny {
protected:
	Node* root;	// the root of the phylogeny, even if the phylogeny is treated as unrooted.
//	std::map<Node*, Node*> anc;	// a map of
	std::vector<Node*> leaves;	// the list of leaf taxa as they appear in the Newick formatted tree
		/**
		 * Each node needs to be stored here in the order in which it appears in the Newick format.
		 * This is so it can be displayed correctly in the compressTraverseWrite method.
		 * However we also need a canonical way in which to index the taxa, so orderedLeaves is a map
		 * from indices 0, 1, ...., (n-1) of the nodes in *alphabetical order*.
		 * Each Node instance stores this index also as the member value idx.
		 * Hence the node idx is calculated after the whole tree is read in.
		 */
	std::set<Node*> orderedLeaves;
	std::set<leafset> splits;
	int labelSpace;	// how much does each vertex label need?
	unsigned int numLeaves;
public:
	Phylogeny();
	virtual ~Phylogeny();

	void addAllDescendants(std::vector<leafset>& splits, Node* v);

	void calculateLeaves();
	void calculateLeaves(Node *v);
	void calculateHeights(Node* v);

	void compressTraverseWrite(std::ostream& os);
	void compressTraverseWrite(std::ostream& os, Node* v);

	int getMaxLabelWidth(Node* v);

	void removeRedundantSplits();
	inline void setRoot(Node* r) { root = r; }
};

} /* namespace parsing */

#endif /* UTILITY_PHYLOGENY_H_ */
