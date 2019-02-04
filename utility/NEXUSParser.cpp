/*
 * NEXUSParser.cpp
 *
 *  Created on: 20 Jul 2016
 *      Author: mac
 */

#include <algorithm>
#include <stdio.h>
#include <sstream>
#include "debugging.h"
#include "../src/project.h"
#include "appexception.h"
#include "node.h"
#include "parser.h"
#include "phylogeny.h"
#include "NEXUSParser.h"

using namespace std;
using namespace flatbush;

extern bool _debugging;

namespace parsing {

NEXUSParser::NEXUSParser(const std::string& fileName, Project* pr)
		: Parser(fileName), proj(pr) {
//	parse();
}

void NEXUSParser::parse() {
	/*
	 * Grammar of NEXUS file:
	 * <nexus> ::= "#NEXUS" { <nexusblock> }
	 * <nexusblock> ::= "begin" <blockname> { <setting> } ( "end" | "endblock" ) ";"
	 * <blockname> ::= <string>
	 * <setting> ::= <id> [ "=" ] <val> [","]
	 */
	TL.reset();
	eat('#');
	eat("nexus");
	while (TL.hasNext()) {
		parseNEXUSBlock();
	}
}

void NEXUSParser::parseBranchLength(Node* v) {
	if (matches(':')) {
		advance();
		double d = getDouble();
		v->setBranchLength(d);
		advance();
	}
}

void NEXUSParser::parseDataBlock() {
try {
	/**
	 * Example:
begin data;
  dimensions ntax=5 nchar=54;
  format datatype=dna missing=? gap=-;
  matrix
    Ephedra       TTAAGCCATGCATGTCTAAGTATGAACTAATTCCAAACGGTGAAACTGCGGATG
    Gnetum        TTAAGCCATGCATGTCTATGTACGAACTAATC-AGAACGGTGAAACTGCGGATG
    Welwitschia   TTAAGCCATGCACGTGTAAGTATGAACTAGTC-GAAACGGTGAAACTGCGGATG
    Ginkgo        TTAAGCCATGCATGTGTAAGTATGAACTCTTTACAGACTGTGAAACTGCGAATG
    Pinus         TTAAGCCATGCATGTCTAAGTATGAACTAATTGCAGACTGTGAAACTGCGGATG
                 [----+--10|----+--20|----+--30|----+--40|----+--50|----]
  ;
end;
	 *
	 * EBNF:
	 * <datablock> ::= "begin" "data" ';' { <component> } <endblock>
	 * <component> ::= ( <dim> | <format> | <matrix> )
	 * <dim>       ::= "dimensions" { ( <ntax> | <nchar> ) } ';'
	 * <format>    ::= "format" { ( <datatype> | <missing> | <gap> ) } ';'
	 * <matrix>    ::= "matrix" { ( <id> <sequence> ) } ;
	 * <endblock>  ::= ( "end" | "endblock" ) ';'
	 */
	bool _debugging(true);
	DEBUG(cout << "Reading DATA block" << endl);
	eat("data");
	ignore(';');	// Be forgiving to missing semicolon here
//	Alignment *A = proj->getAlignment();
	DEBUG(cout << "First token in DATA block is..." << current() << endl);
	while (hasNext() && !matches({ "end", "endblock" })) {
//		advance();
		if (matches("dimensions")) {
			parseDimensionsComponent();
		} else if (matches("format")) {
			parseDataFormat();
		} else if (matches("matrix")) {
			advance();
			DEBUG(cout << "reading matrix" << endl);
			while (hasNext() && !matches(';')) {
				// read sequence name and then sequence; if the same sequence name comes up again, append the sequence.
				string seqID = getString();
				DEBUG(cout << "seqID = " << seqID << endl)
				advance();
				string seq = getString();
				std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);	// have to convert to uppercase. sheesh.
				DEBUG(cout << "seq = " << seqID << endl)
				advance();
				DEBUG(cout << seqID << "\t" << seq << endl);
				map<string, string>* seqs = proj->getRawAlignment();
				if (proj->hasSequence(seqID)) {	// XXX Maybe maybe move this to a method in Alignment class.
					string newSeq = proj->getRawSequence(seqID);
					newSeq += seq;
					(*seqs)[seqID] = newSeq;
				} else {
					(*seqs)[seqID] = seq;
				}
//				DEBUG(cout << "number of sequences in alignment: " << seqs->size() << endl);
			}
			ignore(';');
		}	// end of matrix parameters (alignment)
//		advance();
	}
	proj->checkSequenceLength();
	DEBUG(cout << "finished DATA block" << endl);
	advance();
	ignore(';');
}
catch (app_exception* e) {
	stringstream sstr(e->what());
	sstr << "Exception in NEXUSParser::parseDataBlock()\n";
	throw new app_exception(sstr.str());
}
}

void NEXUSParser::parseDataFormat() {
	bool _debugging = false;
	advance();
	while (hasNext() && !matches(';')) {
		// read the dimensions:
		if (matches("datatype")) {
			advance();
			ignore('=');
			string dt = getString();
//			proj->setDataType(dt);
			DEBUG(cout << "Ignoring datatype field format = " << dt << endl);
			advance();
		} else if (matches("missing")) {
			advance();
			ignore('=');
			proj->setAlignmentMissingChar(getChar());
			advance();
		} else if (matches("gap")) {
			advance();
			ignore('=');
			proj->setAlignmentGapChar(getChar());
			advance();
		}
	}
	ignore(';');	// end of format parameters
	DEBUG(cout << "finished format sub-block" << endl);
}

void NEXUSParser::parseDimensionsComponent() {
	advance();
//	Alignment *A = proj->getAlignment();
	while (hasNext() && !matches(';')) {
		// read the dimensions:
		if (matches("ntax")) {
			advance();
			ignore('=');
			proj->setNTax(getInt());
			advance();
		} else if (matches("nchar")) {
			advance();
			ignore('=');
			proj->setSeqLength(getInt());
			advance();
		}
	}
	ignore(';');
}

enum splitInputFormat {
	split_input_taxonset,
	split_input_binary
};

void NEXUSParser::parseFlatbushBlock() {
	bool _debugging = false;
	eat("flatbush");
	ignore(';');
	splitInputFormat sif(split_input_binary);
	while (hasNext() && !matches({ "end", "endblock" })) {
//		advance();
		if (matches("tree")) {
			Phylogeny* T = new Phylogeny();
			parseNewickFormatTree(T);	// XXX should move this to Phylogeny class.
			T->compressTraverseWrite(cout);
			ignore(';');
		} else if (matches("splits")) {
			advance();

			if (matches({"taxonset", "taxa"}))
				sif = split_input_taxonset;
			string splitStr;
			string s;
			set<string> taxSet;
			switch (sif) {
				case split_input_binary:
					while (!matches(';')) {
						// read a string of 0s and 1s for the split
						advance();
						splitStr = getString();
						proj->addKnownSplit(splitStr);
						advance();
						ignore(',');
					}
					break;
				case split_input_taxonset:
					while (!matches(';')) {
						advance();
						eat('{');
						while (!matches('}')) {
							s = getString();
							DEBUG(cout << s);
							taxSet.insert(s);
							advance();
							DEBUG(cout << "\ttoken = " << current() << endl);
						}
						DEBUG(cout << "reached end of set with token " << current() << endl);
						advance();
						ignore(',');
						DEBUG(cout << "after ignore(,), current token is " << current() << endl);
						DEBUG(cout << "{ "; for (string s : taxSet) cout << s << " "; cout << "}");
						proj->addSplit(taxSet);
						DEBUG(cout << "Added split" << endl);
						DEBUG(cout << "Number of splits stored is " << proj->getSplits()->size());
					}
					break;
			}
//			while (!matches(';')) {
//				DEBUG(cout << "token = " << current() << endl);
//				if (matches({"taxonset", "taxa"})) {
//					// read a set of taxon names; this is half of the split
//					advance();
//					eat('{');
//					string s;
//					set<string> taxSet;
//					while (!matches('}')) {
//						s = getString();
//						DEBUG(cout << s);
//						taxSet.insert(s);
//						advance();
//						DEBUG(cout << "\ttoken = " << current() << endl);
//					}
//					DEBUG(cout << "reached end of set with token " << current() << endl);
//					advance();
//					ignore(',');
//					DEBUG(cout << "after ignore(,), current token is " << current() << endl);
//					DEBUG(cout << "{ "; for (string s : taxSet) cout << s << " "; cout << "}");
////					DEBUG(cout << "taxSet = " << taxSet << endl);
//					proj->addSplit(taxSet);
//					DEBUG(cout << "Added split" << endl);
//					DEBUG(cout << "Number of splits stored is " << proj->getSplits()->size());
//				} else if (matches("binary")) {
//					// read a string of 0s and 1s for the split
//					advance();
//					string splitStr = getString();
//					proj->addKnownSplit(splitStr);
//					advance();
//					ignore(',');
//				}
//			}
		} else if (matches({"outputformat", "outfmt", "-of"})) {
			advance();
			ignore('=');
			if (this->isInt()) {
				int n = getInt();
				switch (n) {
					case fmt_binary_csv:
						proj->outfmt = fmt_binary_csv;
						break;
					case fmt_binary_tabbed:
						proj->outfmt = fmt_binary_tabbed;
						break;
					case fmt_sets:
						proj->outfmt = fmt_sets;
						break;
					default:
						break;
				}
			} else if (this->isString()) {
				if (matches({"binary_csv","binary+csv"})) {
					proj->outfmt = fmt_binary_csv;
				} else if (matches({"binary_tabbed", "binary+tabbed"})) {
					proj->outfmt = fmt_binary_tabbed;
				} else if (matches({"sets", "long"})) {
					proj->outfmt = fmt_sets;
				}
			}
		} else if (matches({"patternfrequencies", "patternfreqs", "pf",
				"patternweights", "pw"})) {
			advance();
			PatternWeights* W = proj->getPatternWeights();
			W->clear();
			string pattern;
			while (!matches(';')) {
				DEBUG(cout << "token = " << current() << endl);
				pattern = getString();
				advance();
				double weight = getDouble();
				W->insert(make_pair(CompactSequence(pattern), weight));
				ignore(',');
				advance();
			}
		}
		if (hasNext()) {
			advance();
		}
	}
	DEBUG(cout << "Completed FLATBUSH block" << endl);
	advance();
	ignore(';');
}

void NEXUSParser::parseMatrix() {
	advance();
//	Alignment *A = proj->getAlignment();
	while (hasNext() && !matches(';')) {
		// read sequence name and then sequence; if the same sequence name comes up again, append the sequence.
		string seqID = getString();
		advance();
		proj->addTaxon(seqID);
		ignore(',');
	}
	ignore(';');
}

void NEXUSParser::parseNEXUSBlock() {
	skipComments();
	eat("begin");
	if (matches("taxa")) {
		parseTaxaBlock();	// now we're linked to a Project, we can add taxon list
	} else if (matches("data")) {
		parseDataBlock();	// put alignment data in here
	} else if (matches("flatbush")) {
		parseFlatbushBlock();
	} else if (matches("sets")) {
		parseSetsBlock();
	} else {
		skipBlock();
	}
}

void NEXUSParser::parseNewickFormatTree(Phylogeny* T) {
	/**
	 * Newick format grammar, adapted from http://evolution.genetics.washington.edu/phylip/newick_doc.html
	 *
	 * <tree>            ::= <descendant_list> [ <rootlabel> ] [ ':' <branchlength> ] ';'
	 * <descendant_list> ::= '(' <subtree> { ',' <subtree> } ')'
	 * <subtree>         ::= (
	 *                          <descendant_list> [ <internal_label> ] [ ':' <branchlength> ]
	 *                       |
	 *                          <leaf_label> [ ':' <branchlength> ]
	 *                       )
	 * <rootlabel>       ::= <label>
	 * <internal_label>  ::= <label>
	 * <leaf_label>      ::= <label>
	 * <label>           ::= <string>
	 * <branchlength>    ::= <number>
	 *
	 * Further, the Newick format can be preceded by an optional tree name:
	 * <treename>        ::= "tree" [ '=' ] <string>
	 */
	if (matches("tree")) {
		advance();
		string treeName = getString();
		advance();
		ignore('=');
	}
	Node* root = new Node();
	parseNewickSubtree(root);
	T->setRoot(root);
	T->calculateLeaves();
	T->calculateHeights(root);
}

void NEXUSParser::parseNewickSubtree(Node* v) {
	/**
	 * See NEXUSParser::parseNewickFormatTree for Newick format
	 */
	DEBUG (cout << current() << endl);

	if (matches('(')) {
		// internal node
		advance();
		Node* child = new Node();
		parseNewickSubtree(child);
		if (isString()) {
			child->setLabel(getString());
			proj->addTaxon(child->getLabel());
			advance();
		} else {
			// create an internal label
		}
		parseBranchLength(child);
		child->setParent(v);
		v->addChild(child);
		while (matches(',')) {
			// siblings!
			advance();
			Node* sib = new Node();
			if (isString()) {
				sib->setLabel(getString());
				proj->addTaxon(sib->getLabel());
				advance();
			} else {
				// create an internal label
			}
			parseNewickSubtree(sib);
			parseBranchLength(sib);
			sib->setParent(v);
			v->addChild(sib);
		}
		eat(')');
	} else if (isString()) {
		Node* leaf = new Node();
		leaf->setLabel(getString());
		proj->addTaxon(leaf->getLabel());
		advance();
		parseBranchLength(leaf);
		leaf->setParent(v);
		v->addChild(leaf);
	}

}

void NEXUSParser::parseSetsBlock() {
	/**
	 * Example:
begin sets;
  CHARSET COI=1-688;
  CHARSET 16S=689-1250;
  CHARSET morph=1251-1322;

  CHARSET COIpos1=1-688\3;
  CHARSET COIpos2=2-688\3;
  CHARSET COIpos3=3-688\3;

  TAXSET outgroup=taxon1 taxon2 taxon3;
  TAXSET NoMorph=taxon33 taxon38 taxon50;
  TAXSET COIonly=1-33;
  TAXSET beetles=22 25 27 33 35 40;
END;
	 */
//	expect("sets");
//	ignore(';');
//	while (matches({"charset", "taxset"})) {
//		if (matches("charset")) {
//			advance();
//
//		}
//	}

}
void NEXUSParser::parseTaxaBlock() {
	/**
	 * Example:
begin taxa;
  dimensions ntax=5;
  taxlabels
    Giardia
    Thermus
    Deinococcus
    Sulfolobus
    Haobacterium
  ;
end;
	 *
	 * EBNF:
	 * <taxablock> ::= "taxa" ';' { <component> } ( "endblock" | "end" ) ';'
	 * <component> ::= <dim> | <taxlabels>
	 * <dim>       ::= "dimensions" { ( <ntax> | <nchar> ) } ';'
	 * <taxlabels> ::= "taxlabels" { <seqID> } ';'
	 *
	 */
	bool _debugging = true;
	DEBUG(cout << "Parsing Taxa block" << endl);
	eat("taxa");
	ignore(';');
//	Alignment *A = proj->getAlignment();
	while (hasNext() && !matches({ "end", "endblock" })) {
//		advance();
		if (matches("dimensions")) {
			parseDimensionsComponent();
		} else if (matches("taxlabels")) {
			while (hasNext()) {
				advance();
				if (matches(';')) {
					break;
				}
				string seqID = getString();
				proj->addTaxon(seqID);
			}
			ignore(';');
		} else if (matches("matrix")) {
			advance();
			while (hasNext() && !matches(';')) {
				// read sequence name and then sequence; if the same sequence name comes up again, append the sequence.
				string seqID = getString();
				advance();
				proj->addTaxon(seqID);
				ignore(',');
			}
			ignore(';');
		} else if (matches("characters")) {
			advance();
			cout << "IMPLEMENT THIS!" << endl;
			if (matches("dimensions")) {
				parseDimensionsComponent();
			} else if (matches("format")) {
				parseDataFormat();
			} else if (matches("matrix")) {
				parseMatrix();
			}
			return;
		}
	}
	proj->gatherTaxonLabels();
	DEBUG(cout << "Number of taxa: " << proj->getNumTaxa() << endl);
	advance();
	ignore(';');

}
void NEXUSParser::skipBlock() {
	while (!matches({"end", "endblock"})) {
		advance();
	}
	advance();
	eat(';');
}

} /* namespace parsing */
