/*
 * project.cpp
 *
 *  Created on: 25 Jul 2016
 *      Author: mac
 */

#include <algorithm>
#include <exception>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <unordered_map>
#include <vector>
#include "project.h"
#include "../utility/appexception.h"
#include "../utility/debugging.h"
#include "../utility/phylogeny.h"
#include "../utility/FASTAParser.h"
#include "../utility/NEXUSParser.h"
#include "../utility/pattern.h"
#include "../utility/CompactSequence.h"
#include "Climb.h"
#include <Eigen/Dense>

using namespace std;

extern bool _debugging;
extern bool _verbose;
extern bool _silent;
extern ofstream out;

/**
 * XXX Put in a data filter:
 * 	remove all sites with any ambiguity codes: can only handle upper case ACGT
 * 	convert all lower case acgt to ACGT respectively
 * 	remove all sites with gaps
 */
//extern const unsigned int k;
// This is a 4x4 representation of the Hadamard matrix with a 1 if the sign is negative and 0 otherwise.
const int negH[4][4] = { { 0, 0, 0, 0 }, { 0, 1, 0, 1 }, { 0, 0, 1, 1 }, { 0, 1, 1, 0 } };

namespace flatbush {

size_t card(const leafset& l) {
	size_t c = 0;
	for (unsigned int i = 0; i < l.size(); ++i) {
		if (l[i]) {
			++c;
		}
	}
	return c;
}

string strmLeafset(const leafset& vec) {
	stringstream ss;
	ss << "{ ";
	for (unsigned int i = 0; i < vec.size(); ++i) {
		if (vec[i]) {
			ss << i << " ";
		}
	}
	ss << "}";
//	for (unsigned int i = 0; i < vec.size(); ++i) {
//		ss << (vec[i]) ? i : '0';
//	}
	return ss.str();
}

string strmLeafsetBinary(const leafset& vec) {
	stringstream ss;
	for (unsigned int i = 0; i < vec.size(); ++i) {
		if (vec[vec.size()-1-i] == true) {
			ss << '1'; //ss << (vec[i]) ? '1' : '0';
		} else {
			ss << '0';
		}
	}
	return ss.str();
}


bool matchesIgnoreCase(string str, const std::vector<std::string>& patterns) {
	// Convert initializer_list to vector and call the (string, vector) overload.
	bool _debugging = false;
	DEBUG(cout << "matchesIgnoreCase(" << str << ", patterns)" << endl);
	transform(str.begin(), str.end(), str.begin(), ::tolower);
	vector<string> vec(patterns.begin(), patterns.end());
//	return matchesIgnoreCase(str, vec);
	for (string s : patterns) {
		transform(s.begin(), s.end(), s.begin(), ::tolower);
		DEBUG(cout << "\tcompared string s=" << s << endl);
		if (str.compare(s) == 0) {
			return true;
		}
	}
	return false;
}

Project::~Project() {}

Project::Project() : parser(nullptr), ft(nullFile), dType(alignmentData),
		ntax(0), seqLength(0), gapChar('-'),
		missingChar('?'), numStates(4), _hasPatternWeights(false),
		outfmt(fmt_binary_csv), _doAllSplits(false) {

}

Project::Project(std::string filename) : parser(nullptr), outfmt(fmt_binary_csv),
		_doAllSplits(false) {
	// load the file with this filename, tokenize it, then parse & create all the defined objects
	read(filename);
}

void Project::addKnownSplit(std::string s) {
	// add a split that is to be checked against, for situations where
	// we know the true tree.
	knownSplits.insert(s);
}

void Project::addRawSequence(string id, string seq) {
	if (seqLength == 0) {
		seqLength = seq.length();
	} else if (seq.length() != seqLength) {
		throw new app_exception("Cannot add sequences of different length to alignment");
	}
	rawAlignment.insert(make_pair(id, seq));
	ntax = rawAlignment.size();
}

void Project::addSequence(string id, string seq) {
	if (seqLength == 0) {
		seqLength = seq.length();
	} else if (seq.length() != seqLength) {
		throw new app_exception("Cannot add sequences of different length to alignment");
	}
	compactAlignment.insert(make_pair(id, seq));
	ntax = rawAlignment.size();
}

/**
 * Add a split in the form of a set of taxon names.
 * @pre the taxon list must be defined! They will be sorted in lex order.
 * The taxon names are accessible from the Alignment, A.
 */
void Project::addSplit(std::set<std::string> S) {
	/**
	 * Add a split by finding the taxon name of each string in the set S and creating a leafset (vector<bool>) for each split.
	 */
	bool _debugging = true;
	leafset split;
	split.resize(ntax);
	if (leafLabels.size() == 0) {
		gatherTaxonLabels();
	}
	DEBUG(cout << "S.size() = " << S.size() << endl);
	DEBUG(cout << "S = { "; for (auto iter = S.begin(); iter != S.end(); ++iter) cout << *iter << " "; cout << "}\n");
	for (set<string>::const_iterator iter = S.begin(); iter != S.end(); ++iter) {
		string str = *iter;
		if (leafLabels.find(str) == leafLabels.end()) {
			stringstream ss;
			for (auto liter = leafLabels.begin(); liter != leafLabels.end(); ++liter) {
				cout << "leaf : " << liter->first << endl;
			}
			ss << "Cannot find leaf named \'" << str << "\' in set of leaves. Sorry, I have to quit.";
			throw new app_exception(ss.str());
		}
		size_t idx = leafLabels.at(str);
		split[idx] = true;
	}
	splits.push_back(split);
}
void Project::addSplit(std::string s) {
	/**
	 * s should be a string of 0s and 1s: convert each 1 at position i to a TRUE in the leafset
	 * Taxa are indexed from right to left, so the FIRST (lex-ordered) taxon name
	 * is the RIGHTMOST / LEAST SIGNIFICANT bit.
	 */
	leafset split;
//	unsigned int ntax = this->getAlignment()->getNumTaxa();
	split.resize(ntax);
	for (size_t i = 0; i < ntax; ++i) {
		switch (s[i]) {
			case '0':
				break;
			case '1':
				split[ntax-1-i] = true;
				break;
			default:
				throw new app_exception("Attempting to add a split in binary string format, but encountered unexpected character.");
		}
	}
	splits.push_back(split);
}

bool Project::adjacent(const leafset& a, const leafset& b) {
	if (card(a) != card(b)) {
		return false;
	}
	return (HammingDistance(a, b) == 1);
}

PatternCounts Project::bootstrapPatternCounts() const {
	PatternCounts sample;
	// not doing anything yet
	return sample;
}

void Project::calcAllSubflatteningErrors() {
//	bool _debugging(true);
//	unsigned int ntax = A.getNumTaxa();
//	unsigned int seqLength = A.getSequenceLength();
//	bool _verbose = true;
//	bool _debugging = true;
	leafset right;
	vector<Pattern<char>> rows;
	vector<Pattern<char>> columns;
	Pattern<char> pat(ntax, '.');
	Pattern<unsigned int> idx;
	Pattern<char> nullPattern(ntax, ' ');
	if (splits.size() == 0) {
		if (_verbose) {
			cout << "No splits to test! Better do them all." << endl;
		}
		gatherAllNontrivialHalfSplits(_verbose);
	} else if (_doAllSplits) {
		if (_verbose) {
			cout << "--do-all-splits command-line flag invoked so doing the lot." << endl;
		}
		gatherAllNontrivialHalfSplits(_verbose);
	}
	if (_verbose) {
		cout << "Known splits:\n";
		for (auto iter = knownSplits.begin(); iter != knownSplits.end(); ++iter)
			cout << "\t" << *iter << endl;
		cout << endl;
	}
	if (out.good()) {
		switch (outfmt) {
			case fmt_binary_csv:
				out << "error, left, right, intree, size" << endl;
				break;
			case fmt_binary_tabbed:
				out << "error\tleft\tright\tintree\tsize" << endl;
				break;
			case fmt_sets:
				out << "error, left, right, intree, size" << endl;
				break;
			default: break;
		}
	}
	for (leafset left : splits) {
//		DEBUG(cout << "left = " << strmLeafset(left) << '\t');
		right = left;
		right.flip();
//		DEBUG(cout << "right = " << strmLeafset(right) << endl);
		putHalfSplitMask(rows, left);
		putHalfSplitMask(columns, right);

		Eigen::MatrixXd F(rows.size(), columns.size());
//		Eigen::SparseMatrix<double> SM(rows)
		double error = 0.0;
		if (_verbose) {
			cout << endl << "Subflattening for split " << setw(ntax) << strmLeafset(left)
					<< ":" << strmLeafset(right) << ", " << describeSplit(left) << endl;
		}
		unsigned int r = 0;
		unsigned int c = 0;
//		cout << setw(static_cast<int>(log(A.getSequenceLength()))+15);
		if (_verbose) {
			cout << "\t" << nullPattern << "\t";
			for (auto col : columns) {
				cout << setw(10-ntax) << col << "\t";
			}
			cout << endl;
		}
		for (auto row : rows) {
			c = 0;
			if (_verbose) {
				cout << "\t" << row << "\t";
			}
			for (auto col : columns) {
				combine(pat, row, col);
				encode(idx, pat);
				F(r, c) = countSignedSum(idx);
				F(r, c) /= getSequenceLength();
				if (_verbose) {
					cout << setw(11) << F(r,c) << "\t";
				}
				++c;
			}
			if (_verbose) {
				cout << endl;
			}
			++r;
		}
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(F);
		Eigen::MatrixXd svalues = svd.singularValues().transpose();
		if (_verbose) {
			cout << "\tSingular Values found by Jacobi SVD: " << svalues << endl;
		}
		/**
		 * ERROR calculation:
		 * Sum from k to the number of columns of F the square of the singular value
		 */
		error = 0.0;
		for (unsigned int erIdx = numStates; erIdx < svalues.cols(); ++erIdx) {
			error += svalues(erIdx) * svalues(erIdx);
		}
		double normer = error;
		for (unsigned int i = 0; i < numStates; ++i) {
			normer += svalues(i) * svalues(i);
		}
		error = sqrt(error / normer);
		if (_verbose) {
			cout << "\tError from rank " << numStates << " : " << error << endl;
		}
		if (out.good()) {
//			out.precision(6);
//			out.setf(std::ios::fixed, std:: ios::floatfield);
			// TODO In this part, output a flag too, to say whether the split is in the tree.
//			DEBUG(cout << idx << "\t" << error << endl);
			string leftString = strmLeafsetBinary(left);
			string rightString = strmLeafsetBinary(right);
			DEBUG(cout << setw(8) << error << ", " << leftString << ", " << rightString << endl);
//			cout << "split: " << leftString << "\t" << rightString << endl;
			int _inTree = 0;
			if (knownSplits.find(leftString) != knownSplits.end()
					|| knownSplits.find(rightString) != knownSplits.end()) {
				_inTree = 1;
			};
			unsigned int size = card(left);
			if (2*size > ntax) {
				size = ntax - size;
			}
			switch (outfmt) {
				case fmt_binary_csv:
					out << setw(8) << error << ", " << leftString << ", "
							<< rightString << ", " << _inTree << ", " << size << endl;
					break;
				case fmt_binary_tabbed:
					out << setw(8) << error << '\t' << leftString << '\t'
							<< rightString << '\t' << _inTree << '\t' << size << endl;
					break;
				case fmt_sets:
					out << setw(8) << error << ", " << strmLeafset(left)
							<< ", " << strmLeafset(right) << _inTree << ", " << size << endl;
					break;
				default: break;
			}
		}
//		if (_verbose) {
//		cout << error << ", " << strmLeafset(left) << ", " << strmLeafset(right) << endl;
		if (!_silent) {
			cout << error << ", " << strmLeafset(left) << ":" << strmLeafset(right) << ", " << describeSplit(left) << endl;
		}
//		}

	}

}

/**
 * Calculate the pattern weights if they're not available.
 */
void Project::calcPatternWeights() {
	if (_hasPatternWeights) {
		return;
	}
	/**
	 * Here I need to choose between using the alignment and using a list of splits with weights.
	 */
	switch (getDataType()) {
	case alignmentData:
		compressAlignment();
		countPatterns(); // puts counts into patternCounts
		break;
	case patternWeightData:
		// should check the pattern weights are, e.g., all non-negative?
		try {
			for (auto iter = patternWeights.begin(); iter != patternWeights.end(); ++iter) {
				if (iter->first.getLength() != ntax) {
					throw new runtime_error("Site pattern length disagrees with ntax");
				}
				if (iter->second < 0.0) {
					stringstream ss;
					ss << "Negative weight for pattern " << iter->first << ", which is prohibited.";
					throw new runtime_error(ss.str());
				}
			}
		} catch (runtime_error& re) {
			cerr << "Runtime Error in Project::calcPatternWeights: " << re.what() << endl;
		}
		break;
	default:
		break;
	};
	_hasPatternWeights = true;
}

void Project::checkSequenceLength(bool nothrow) {
	if (compactAlignment.size() == 0) {
		return;	// nothing to check as there are no sequences
	}
	for (auto iter : compactAlignment) {

	}
	auto iter = compactAlignment.begin();
	unsigned int sl = iter->second.getLength();
	++iter;
	while (iter != compactAlignment.end()) {
		if (iter->second.getLength() != sl) {
			stringstream msg("two sequences of different length: sequence \'");
			msg << iter->first << "\' has length " << iter->second.getLength() << " but " << sl << " was expected.";
			throw new app_exception(msg.str());
		}
		++iter;
	}
	if (nothrow) {
		seqLength = sl;
		return;
	}
	if (sl != seqLength) {
		stringstream msg("seqlength was found to be ");
		msg << sl << " but " << seqLength << " was expected.";
		throw new app_exception(msg.str());
	}
}

void Project::compressAlignment() {
	/**
	 * Search each character in each sequence
	 * If a character is undesirable, flag that site for removal
	 */
	vector<bool> _uselessSite(seqLength, false);
	string acceptable("ACGT");
	unsigned int numRemoved(0);
	for (auto iter = rawAlignment.begin(); iter != rawAlignment.end(); ++iter) {
		string& seq = iter->second;
		for (unsigned int i = 0; i < seqLength; ++i) {
			if (acceptable.find(seq[i]) == string::npos) {
//				cout << "cannot find character " << seq[i] << " in " << acceptable << endl;
				_uselessSite[i] = true;
			}
		}
	}
	DEBUG(
		for (unsigned int i = 0; i < seqLength; ++i) {
			if (_uselessSite[i]) {
				cout << "Removing this site:\t";
			} else {
				cout << "keeping this site:\t";
			}
			for (auto iter = rawAlignment.begin(); iter != rawAlignment.end(); ++iter) {
				cout << iter->second[i];
			}
			cout << endl;
		}
	);
	for (auto iter = rawAlignment.begin(); iter != rawAlignment.end(); ++iter) {
		stringstream seq;
		for (unsigned int i = 0; i < seqLength; ++i) {
			if (!_uselessSite[i]) {
				seq << iter->second[i];
			}
		}
		iter->second = seq.str();
	}
	// count the number of useless sites:
	for (unsigned int i = 0; i < seqLength; ++i) {
		numRemoved += (_uselessSite[i]) ? 1 : 0;
	}
	if (numRemoved > 0) {
		cout << "Removed " << numRemoved << " sites with gaps or other uselessness." << endl;
		cout << "New sequence length = " << seqLength << endl;
	}
	/**
	 * Now convert the strings into CompactSequence objects:
	 */
	compactAlignment.clear();
	for (auto iter : rawAlignment) {
		compactAlignment.insert(pair<string, CompactSequence>(iter.first, CompactSequence(iter.second)));
	}
}

/**
 * Search the Alignment A to find all the names of all the sequences;
 * then put them in a map (in order);
 * then assign each an index, which will be used in pattern counting (and possibly elsewhere).
 */
void combine(Pattern<char>& result, const Pattern<char> x, const Pattern<char> y) {
	result.resize(x.size());
	unsigned int numDots = 0;
	for (unsigned int i = 0; i < x.size(); ++i) {
		numDots = (x[i] == '.') ? 1 : 0;
		numDots += (y[i] == '.') ? 1 : 0;
		if (numDots != 1) {
			throw new app_exception("Cannot combine these patterns: must have exactly one \'.\' at the same position");
		}
		result[i] = (x[i] == '.') ? y[i] : x[i];
	}
}

void Project::countAllSignedSums(bool _verbose) {
	/**
	 * Count all the signed sums that are possible.
	 * The possible patterns for signed sums are constrained as follows:
	 *
	 * 1. The patterns are of length n
	 * 2. Each pattern has at most 2 non-0 entries
	 *
	 * This function iterates over all possible such patterns and accumulates their signed counts in one map,
	 * sequencePatternCount.
	 */
	if (_verbose) {
		cout << "countAllSignedSums()..." << endl;
	}
	unsigned int k = 4;
	int signedSum = 0;
	Pattern<unsigned int> pat;
	int numPatterns = 0;
	// all 0s:
	for (unsigned int i = 0; i < ntax; ++i) {
		pat.append(0);
	}
	++numPatterns;
	signedSum = countSignedSum(pat);
	signedCountSums.insert(make_pair(pat, signedSum));
	if (_verbose) {
		cout << pat << " : " << signedSum << endl;
	}
	// singletons:
	for (unsigned int i = 0; i < ntax; ++i) {
		for (unsigned int a = 1; a < k; ++a) {
			// fill in the pattern:
			pat = 0;
			pat[i] = static_cast<unsigned int>(a);
			++numPatterns;
			signedSum = countSignedSum(pat);
			signedCountSums.insert(make_pair(pat, signedSum));
			if (_verbose) {
				cout << pat << " : " << signedSum << endl;
			}
		}
	}
	// matched doubles:
	for (unsigned int i = 0; i < ntax-1; ++i) {
		for (unsigned int j = i+1; j < ntax; ++j) {
			for (unsigned int a = 1; a < k; ++a) {
				pat = 0;
				pat[i] = a;
				pat[j] = a;
				++numPatterns;
				signedSum = countSignedSum(pat);
				signedCountSums.insert(make_pair(pat, signedSum));
				if (_verbose) {
					cout << pat << " : " << signedSum << endl;
				}
			}
		}
	}
	// unmatched doubles:
	for (unsigned int i = 0; i < ntax-1; ++i) {
		for (unsigned int j = i+1; j < ntax; ++j) {
			for (unsigned int a = 1; a < k-1; ++a) {
				for (unsigned int b = a+1; b < k; ++b) {
					pat = 0;
					pat[i] = a;
					pat[j] = b;
					++numPatterns;
					signedSum = countSignedSum(pat);
					signedCountSums.insert(make_pair(pat, signedSum));
					if (_verbose) {
						cout << pat << " : " << signedSum << endl;
					}
					pat[j] = a;
					pat[i] = b;
					++numPatterns;
					signedSum = countSignedSum(pat);
					signedCountSums.insert(make_pair(pat, signedSum));
					if (_verbose) {
						cout << pat << " : " << signedSum << endl;
					}
				}
			}
		}
	}
	if (_verbose) {
		cout << "Total number of patterns = " << numPatterns << endl;
	}
}

void Project::countPatterns() {
	bool _debugging(false);
	DEBUG(cout << "Alignment::countPatterns()" << endl);
//	ntax = compactAlignment.size();
	DEBUG(cout << "numTaxa = " << ntax << endl);
	DEBUG(cout << "size of rawAlignment = " << rawAlignment.size() << endl);
	patternCounts.clear();
	DEBUG(cout << "bugfinding... made it to here.\n")
	char *pat = new char[ntax+1];
	for (unsigned int i = 0; i < seqLength; ++i) {
		unsigned int j(0);
		for (auto &t : rawAlignment) {
			pat[j] = t.second[i];
			++j;
		}
		pat[j] = '\0';
		string str(pat);
		DEBUG(cout << "str = " << str << endl);
		patternCounts[str]++; // doing an implicit type conversion here.. worth a try.
	}
	delete [] pat;
}


int Project::countSignedSum(Pattern<unsigned int> pat) {
/**
 * @para patFreq is the pattern frequencies as unsigned integers; \f$ \mathcal{C} \f$ in the accompanying paper.
 * sigma = 1 is our "magic" state
 * pat is the (m+n)-tuple of numbers in the range [1,k], with at most two entries that are not equal to k:
 * E.g., with n=4, k=3, these are permitted:
 * 	3333 3331 3313 3133 1333 3332 3323 3233 2333 3311 3131 3113 ...
 * and these are not:
 * 	3111 2222
 * The transformed count is:
 \f$
  	F'_{x_1 x_2 x_3 ...x_{m+n}} = \sum_{j_1,...,j_{m+n}\in \Sigma} H_{x_1 j_1} H_{x_2 j_2} ... H_{x_{m+n} j_{m+n}} \mathcal{C}_{j_1 j_2 ... j_{m+n}}
 \f$
 * There are
 */
	// memoisation:
	if (signedCountSums.find(pat) != signedCountSums.end()) {
//		cout << "pattern already found -- the signed sum is " << (*signedSums.find(pat)).second << endl;
		return (*signedCountSums.find(pat)).second;
	}
	int count = 0;
	unsigned int patLength = ntax;
	for (auto iter = patternCounts.begin(); iter != patternCounts.end(); ++iter) {
		CompactSequence p = iter->first;
//		cout << p << ", ";
		int term = 0;
		for (unsigned int i = 0; i < patLength; ++i) {
			if (pat[i] > 3) {
				throw new app_exception("pat element is out of range");
			}
			if (p[i] != 'C' && p[i] != 'A' && p[i] != 'T' && p[i] != 'G') {
				stringstream sstr;
				sstr << "p = " << p << "; p[" << i << "] = " << p[i] << " is not in {A,C,G,T}" << endl;
				throw new app_exception(sstr.str());
			}
			if (encoding.at(p[i]) < 0 || encoding.at(p[i]) > 3) {
				throw new app_exception("encoding is out of range [0,3]");
			}
			term += negH[pat[i]][encoding.at(p[i])];
//			cout << term;
		}
		term = (term % 2 == 0) ? 1 : -1;	// If number of negative 1s in the product is even, then the factor becomes 1; else it's -1.
//		cout << term << " * #" << p << " = " << term * static_cast<int>((*iter).second) << endl;
		term *= static_cast<int>(iter->second);
//		if (term == -1) cout << "-";
//		cout << (*iter).second << endl;
		count += term;
	}
	return count;
}

string Project::describeSplit(leafset pat) {
	string in;
	string out;
	int taxon = 0;
	for (auto it : compactAlignment) {
		if (pat[taxon]) {
			in += it.first;
			in += " ";
		} else {
			out += it.first;
			out += " ";
		}
		++taxon;
	}
	string str = "{ ";
	str += in;
	str += "| ";
	str += out;
	str += "}";
	return str;
}

void Project::doLandscapeAnalysis() {
	bool _debugging(true);
	leafset right;
	vector<Pattern<char>> rows;
	vector<Pattern<char>> columns;
	Pattern<char> pat(ntax, '.');
	Pattern<unsigned int> idx;
	Pattern<char> nullPattern(ntax, ' ');
	walks.clear();
	if (splits.size() == 0 || _doAllSplits) {
		gatherAllNontrivialHalfSplits(_verbose);
	}
	calcPatternWeights();
	unsigned int walkLength;
	double score(0.0);
	for (leafset left : splits) {
		score = getSplitError(left);
		leafset finish;
//		cout << "Start: " << strmLeafset(left) << ',' << score << endl;
		if (walks.find(left) == walks.end()) {
//			DEBUG(cout << "New starting point in base loop: " << strmLeafset(left) << endl);
//			steepestDescentFromHere(left, walkLength);
			recursiveSteepestDescent(left, score, finish, walkLength);
		}
//		cout << strmLeafset(left) << ',' << score << ',' << walkLength << endl;
	}
	if (_verbose) {
		cout << "Landscape analysis complete." << endl;
	}
}

void Project::exportLandscape(string & filename) {
	/**
	 * Create a dot file (for GraphViz) with each node being a split,
	 * 	edges between splits if they're adjacent under the current perturbation,
	 * 	and a weight for each split corresponding to its score.
	 * Not sure how to convey weights in dot... except by colour maybe?
	 */
//	ofstream out;
//	try {
//		out.open(filename.c_str());
//		if (!out.good()) {
//			return;
//		}
//		out << "graph SplitSpace {\n";
//		for (const leafset& a : splits) {
//			out << "\ts" << strmLeafsetBinary(a) << " [label=\"" << strmLeafset(a) << "\"];";
//		}
//		for (vector<leafset>::const_iterator i = splits.begin(); i != splits.end(); ++i) {
//			for (vector<leafset>::const_iterator j = i+1; j!= splits.end(); ++j) {
//				if (adjacent(*i, *j)) {
//					cout << "\ts" << strmLeafsetBinary(*i) << " -- s" << strmLeafsetBinary(*j) << ";";
//				}
//			}
//		}
//		out << "}\n";
//	}
//	catch (exception& e) {
//		cerr << "exporting landscape failed." << e.what() << endl;
//		return;
//	}
}

CompactSequence Project::getSequence(const string& seqID) {
	return compactAlignment.at(seqID);
}

bool Project::hasSequence(string seqID) {
	return (compactAlignment.find(seqID) != compactAlignment.end());
}

//void Project::setDataType(const string& dt) {
//	dType = dt;
//	if (!strcmp(dt.c_str(), "DNA")) {
//		numStates = 4;
//	} else if (!strcmp(dt.c_str(), "AA")) {
//		numStates = 20;
//	} else if (!strcmp(dt.c_str(), "codon")) {
//		numStates = 64;
//	} else if (!strcmp(dt.c_str(), "binary")) {
//		numStates = 2;
//	} else {
//		throw new app_exception("Cannot understand this data type: please choose DNA, AA, codon, or binary (case sensitive)");
//	}
//	if (numStates != 4) {
//		throw new app_exception("Can only handle four states at present: please choose DNA.");
//	}
////	if (parsing::matchesIgnoreCase(dataType, {"DNA", "RNA"})) {
////		numStates = 4;
////	} else if (parsing::matchesIgnoreCase(dataType, {"AA", "aminoacid", "protein"})) {
////		numStates = 20;
////	} else if (parsing::matchesIgnoreCase(dataType, {"codon", "codons"})) {
////		numStates = 64;
////	} else {
////		stringstream msg;
////		msg << "Alignment::setDataType(" << dataType << "). Cannot identify the data type.";
////		throw new app_exception(msg.str());
////	}
//}

void Project::showPatternCounts() {
	for (auto &p : patternCounts) {
		cout << p.first << "\t" << p.second << endl;
	}
	cout << "Total number of unique patterns found: " << patternCounts.size() << endl;
}

void Project::showSignedCounts() {
	for (pair<Pattern<unsigned int>, int> p : signedCountSums) {
		cout << p.first << "\t" << p.second << endl;
	}
}

// Alphabetical order of character states A C G T by default:
//but just out of interest....
const std::map<char, int> Project::encoding = {
		std::make_pair('A', 0),
		std::make_pair('C', 1),
		std::make_pair('G', 2),
		std::make_pair('T', 3),
		std::make_pair('0', 0),
		std::make_pair('1', 1),
		std::make_pair('2', 2),
		std::make_pair('3', 3)
	};

void Project::gatherAllNontrivialHalfSplits(bool _verbose) {
//	bool _debugging = true;
//	int ntax = A.getNumTaxa();
	if (ntax > 20) {
		throw new app_exception("that's too many splits, sorry. Stick to 20 taxa or fewer for this method.");
	}
	splits.clear();
	leafset half(ntax, false);
	for (int i = 3; i < (1<<(ntax-1)); ++i) { // XXX At some point a Gray code would be cool to implement for this.
		// check out coding answer here: https://stackoverflow.com/questions/17490431/gray-code-increment-function
		int card = 0;
		for (int t = 0; t < ntax-1; ++t) { // don't include the last taxon, labelled (n-1).
			if (i & (1<<t)) {
				half[t] = true;
				++card;
			} else {
				half[t] = false;
			}
		}
		if (card <= 1 || card >= ntax-1) {
			continue;
		}
		DEBUG(cout << "adding partition { ";
			for (int x = 0; x < ntax; ++x) {
				if (half[x]) {
					cout << x << " ";
				}
			}
			cout << "}" << endl
		);
		splits.push_back(half);
//		if (_verbose) {
//			cout << strmLeafset(half) << " ";
//		}
	}
	cout.flush();
}

void Project::gatherTaxonLabels() {
	for (auto taxSeq : compactAlignment) {
		leafLabels[taxSeq.first] = 0;
	}
	int i = 0;
	for (map<string, int>::iterator iter = leafLabels.begin(); iter != leafLabels.end(); ++iter) {
		leafLabels[iter->first] = i;
		++i;
	}
	ntax = leafLabels.size();
	DEBUG(cout << "{ ");
	DEBUG(for (map<string, int>::const_iterator citer = leafLabels.begin(); citer != leafLabels.end(); ++citer) cout << citer->first << " ");
	DEBUG(cout << "}\n");
}

void Project::putAllAdjacentSplits(vector<leafset>& N, const leafset L) {
	bool _debugging(false);
	N.clear();
//	// First perturbation: move a single taxon from one side to the other:
//	for (unsigned int i = 0; i < ntax; ++i) {
//		leafset a(L);
//		a[i] = 1 - a[i]; // flip the ith bit
//		if (isTrivialSplit(a)) {
//			continue;
//		}
//		DEBUG(cout << "original split: " << strmLeafset(L) << "; next split by single flip = " << strmLeafset(a) << endl);
//		N.push_back(a);
//	}
	// Second perturbation: swap two taxa:
	for (unsigned int i = 0; i < ntax; ++i) {
		if (L[i] == 0) {
			continue;
		}
		for (unsigned int j = i+1; j < ntax; ++j) {
			leafset a(L);
			a[i] = 0;
			if (a[j]) {
				continue;
			}
			a[j] = 1;
			DEBUG(cout << "original split: " << strmLeafset(L) << "; next split by double flip = " << strmLeafset(a) << endl);
			if (isTrivialSplit(a)) {
				continue;
			}
			if (a[ntax-1]) {
				// replace with the complement as this set includes (ntax-1)
				for (unsigned int pos = 0; pos < a.size(); ++pos) {
					a[pos] = 1 - a[pos];
				}
			}
			N.push_back(a);
		}
	}
}

leafset Project::getRandomSplit() {
	leafset L(ntax);
	do {
		for (unsigned int i = 0; i < ntax; ++i) {
			L[i] = (drand48() < 0.5) ? 1 : 0;
		}
	} while (isTrivialSplit(L));
	return L;
}

double Project::getSplitError(const leafset& left) {
	bool _debugging(false);
	if (splitError.find(left) != splitError.end()) {
		return splitError.at(left);
	}
//	bool _verbose(true);
//	DEBUG(cout << "input split (left) = " << strmLeafset(left) << '\t');
	leafset right(left);
	right.flip();
//	DEBUG(cout << "right = " << strmLeafset(right) << endl);
	vector<Pattern<char>> rows;
	vector<Pattern<char>> columns;
	Pattern<char> pat(ntax, '.');
	Pattern<unsigned int> idx;
	putHalfSplitMask(rows, left);
	putHalfSplitMask(columns, right);

	Eigen::MatrixXd F(rows.size(), columns.size());
//		Eigen::SparseMatrix<double> SM(rows)
	double error = 0.0;
	if (_verbose) {
		cout << endl << "Subflattening for split " << setw(ntax) << strmLeafset(left)
				<< ":" << strmLeafset(right) << ", " << describeSplit(left) << endl;
	}
	unsigned int r = 0;
	unsigned int c = 0;
//		cout << setw(static_cast<int>(log(A.getSequenceLength()))+15);
	if (_verbose) {
		Pattern<char> nullPattern(ntax, ' ');
		cout << "\t" << nullPattern << "\t";
		for (auto col : columns) {
			cout << setw(10-ntax) << col << "\t";
		}
		cout << endl;
	}
	for (auto row : rows) {
		c = 0;
		if (_verbose) {
			cout << "\t" << row << "\t";
		}
		for (auto col : columns) {
			combine(pat, row, col);
			encode(idx, pat);
			F(r, c) = countSignedSum(idx);
			F(r, c) /= getSequenceLength();
			if (_verbose) {
				cout << setw(11) << F(r,c) << "\t";
			}
			++c;
		}
		if (_verbose) {
			cout << endl;
		}
		++r;
	}
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(F);
	Eigen::MatrixXd svalues = svd.singularValues().transpose();
	if (_verbose) {
		cout << "\tSingular Values found by Jacobi SVD: " << svalues << endl;
	}
	/**
	 * ERROR calculation:
	 * Sum from k to the number of columns of F the square of the singular value
	 */
	error = 0.0;
	for (unsigned int erIdx = numStates; erIdx < svalues.cols(); ++erIdx) {
		error += svalues(erIdx) * svalues(erIdx);
	}
	double normer = error;
	for (unsigned int i = 0; i < numStates; ++i) {
		normer += svalues(i) * svalues(i);
	}
	error = sqrt(error / normer);
	if (_verbose) {
		cout << "\tError from rank " << numStates << " : " << error << endl;
	}
	if (out.good()) {
//			out.precision(6);
//			out.setf(std::ios::fixed, std:: ios::floatfield);
		// TODO In this part, output a flag too, to say whether the split is in the tree.
//			DEBUG(cout << idx << "\t" << error << endl);
		string leftString = strmLeafsetBinary(left);
		string rightString = strmLeafsetBinary(right);
		DEBUG(cout << setw(8) << error << ", " << leftString << ", " << rightString << endl);
//			cout << "split: " << leftString << "\t" << rightString << endl;
		int _inTree = 0;
		if (knownSplits.find(leftString) != knownSplits.end()
				|| knownSplits.find(rightString) != knownSplits.end()) {
			_inTree = 1;
		};
		unsigned int size = card(left);
		if (2*size > ntax) {
			size = ntax - size;
		}
		switch (outfmt) {
			case fmt_binary_csv:
				out << setw(8) << error << ", " << leftString << ", "
						<< rightString << ", " << _inTree << ", " << size << endl;
				break;
			case fmt_binary_tabbed:
				out << setw(8) << error << '\t' << leftString << '\t'
						<< rightString << '\t' << _inTree << '\t' << size << endl;
				break;
			case fmt_sets:
				out << setw(8) << error << ", " << strmLeafset(left)
						<< ", " << strmLeafset(right) << _inTree << ", " << size << endl;
				break;
			default: break;
		}
	}
	splitError[left] = error;
	return error;
}

vector<leafset>* Project::getSplits() {
	return &splits;
}

unsigned int Project::HammingDistance(const leafset& a, const leafset& b) {
	/**
	 * Count the number of differences between a and b.
	 */
	unsigned int h(0);
	for (unsigned int i = 0; i < a.size(); ++i) {
		h += (a[i] != b[i]) ? 1 : 0 ;
	}
	return h;
}

bool Project::isTrivialSplit(const leafset& split) {
	unsigned int size = card(split);
//	bool _debugging(true);
	bool _trivial((size < 2) || (size > ntax-2));
//	DEBUG(cout << "testing for trivial split: split = " << strmLeafset(split)
//			<< "; size = " << size << "; trivial = " << _trivial << endl);
	return _trivial;
}

double Project::findBestNeighbouringSplit(leafset& bestNeighbour, leafset L, double score) {
	bool _debugging(false);
	vector<leafset> N;
	DEBUG(cout << "finding all adjacent splits... "; cout.flush());
	putAllAdjacentSplits(N, L);
	DEBUG(cout << N.size() << " splits found.\n");
	double bestScore(score);
	bestNeighbour = L;
	for (auto s : N) {
		double d = getSplitError(s);
		DEBUG(cout << "adjacent split : " << strmLeafset(s) << ", score = " << d << endl);
		if (d < bestScore) {
			bestScore = d;
			bestNeighbour = s;
		}
	}
	DEBUG(cout << "*BEST* split found: " << strmLeafset(bestNeighbour) << ", score = " << bestScore << endl);
	return bestScore;
}

double Project::recursiveSteepestDescent(const leafset& currentSplit, double currentScore, leafset& finish, unsigned int& walkLength) {
	/**
	 * Check all neighbours for that split with the minimal score
	 * Select that split, b
	 * If walk(b) does not exist then
	 * 	recursiveSteepestDescent(b)
	 * walk(current) to contain walks(b).finish(), walks(b).finalScore(), and walks(b).walkLength+1
	 */
	bool _debugging(true);
	leafset b;
//	DEBUG(cout << "recursiveSteepestDescent: " << strmLeafset(currentSplit) << endl);
	double bestScore = findBestNeighbouringSplit(b, currentSplit, getSplitError(currentSplit));
//	DEBUG(cout << "bestScore = " << bestScore << endl);
	if (walks.find(b) != walks.end()) {
		// we have this descent / climb already from b to a local opt, so just tack this on to the currentSplit climb:
		Climb & desc = walks.at(b);
//		DEBUG(cout << "This is a known walk from " << strmLeafset(b) << " to " << strmLeafset(desc.finish)
//				<< "; final score = " << desc.finalScore << " and length " << desc.walkLength << ".\n");
		Climb climb(currentSplit, getSplitError(currentSplit), desc.finish, desc.walkLength+1, desc.finalScore);
//		DEBUG(cout << "New climb: " << strmLeafset(currentSplit) << " -> " << strmLeafset(climb.finish)
//				<< ", score = " << climb.finalScore << " in " << climb.walkLength << " steps.\n");
		walks[currentSplit] = climb;
		walkLength = desc.walkLength+1;
		finish = desc.finish;
		cout << climb << endl;
//		DEBUG(cout << "size of walks map = " << walks.size() << endl);
		return desc.finalScore;
	}
	// Else, we don't have this walk recorded already.
//	DEBUG(cout << "Have not encountered this walk before.\n");
//	DEBUG(cout << "bestScore = " << bestScore << "; currentScore = " << currentScore
//			<< "; difference = " << (bestScore - currentScore) << endl);
	if (b == currentSplit) {
		// If there are no more improving neighbours, this is a walk of length 0
		Climb climb(currentSplit, getSplitError(currentSplit), currentSplit, 0, bestScore);
		finish = currentSplit;
		walkLength = 0;
		walks[currentSplit] = climb;
		cout << climb << endl;
//		DEBUG(cout << "size of walks map = " << walks.size() << endl);
		return bestScore;
	}
	// Otherwise, we have an improving neighbour, so go to it and repeat:
//	if (b == current) {
//		return bestScore;
//	}
//	DEBUG(cout << "Moving on to best neighbour: " << strmLeafset(b) << endl);
	bestScore = recursiveSteepestDescent(b, currentScore, finish, walkLength);
	Climb climb(currentSplit, getSplitError(currentSplit), finish, walkLength+1, bestScore);
	cout << climb << endl;
	++walkLength;
	walks[currentSplit] = climb;
//	DEBUG(cout << "size of walks map = " << walks.size() << endl);
	return bestScore;
}

double Project::steepestDescentFromHere(leafset& current, unsigned int& walkLength) {
	if (walks.find(current) != walks.end()) {
		Climb& climb = walks.at(current);
		walkLength = climb.walkLength + 1; // XXX here
	}
	double score = getSplitError(current);
	DEBUG(cout << "initial score = " << score << endl);
	DEBUG(cout << "#Steepest Descent to find locally optimal splits\n");
	DEBUG(cout << "split, score\n");
	double nextScore(score);
	leafset next(current);
	vector<leafset> neighbours;
	walkLength = 0;
	while (true) {
		DEBUG(cout << strmLeafset(current) << ',' << score << endl);
		nextScore = findBestNeighbouringSplit(next, current, score);
		if (nextScore == score) {
			break;
		}
		++walkLength;
		score = nextScore;
		current = next;
	};
	return score;
}

double Project::steepestDescentRandomStart() {
	bool _debugging(true);
	leafset current = getRandomSplit();
	unsigned int ignoreWalkLength;
	return steepestDescentFromHere(current, ignoreWalkLength);
}

void Project::putHalfSplitMask(vector<Pattern<char>>& idx, flatbush::leafset X) {

	/**
	 * put the "half" split patterns that can go into the row or column indices.
	 * The format is like this:
	 * 	..0..100..
	 * which corresponds to the split (counting from 0, not 1) of { 2, 5, 6, 7 } | { 0, 1, 3, 4, 8, 9 }.
	 * The pattern is the bit that changes:
	 * 	..0..000..
	 * 	..0..001..
	 * 	..0..010..
	 * 	etc.
	 */
	idx.clear();
//	int k = A.getNumStates();
	Pattern<char> pat;

	// all 0s:
	unsigned int n = X.size();
	for (unsigned int i = 0; i < n; ++i) {
		char ch = '.';
		if (X[i]) {
			ch = '0';
		}
		pat.append(ch);
	}
//	cout << pat << endl;
	idx.push_back(pat);
//	unsigned int card = std::count_if(X.begin(), X.end(), bind1st(equal_to<bool>(), true));
	unsigned int card = 0;
	for (unsigned int i = 0; i < n; ++i) {
		if (X[i]) {
			++card;//std::count(X.begin(), X.end(), static_cast<bool>(true));
		}
	}
	vector<unsigned int> position;
	// in the example above, positions will be set to be the vector ( 2, 5, 6, 7 )
	for (unsigned int i = 0; i < n; ++i) {
		if (X[i]) {
			position.push_back(i);
//			cout << "position += " << i << endl;
		}
	}

	// singletons:
	for (unsigned int pos = card-1; ; --pos) {
		pat = '.';
//		cout << pat << endl;
		for (unsigned int a = 1; a < numStates; ++a) {
			for (unsigned int left = 0; left < card; ++left) {
				pat[position[left]] = '0';
			}
			pat[position[pos]] = static_cast<char>('0' + a);
//			cout << pat << endl;
			idx.push_back(pat);
		}
		if (pos == 0) {
			break;
		}
	}
}


void Project::read(const std::string& filename) {
	bool _debugging = true;
	DEBUG(cout << "input filename = " << filename << endl);
//	DEBUG(cout << "finding last position of . in filename..." << endl);
	size_t pos = filename.find_last_of(".");
//	DEBUG(cout << "position is " << pos << endl);
	string suffix = filename.substr(pos);
//	DEBUG(cout << "suffix is " << suffix << endl);
//	DEBUG(cout << "FASTA accepted suffixes: ";
//		vector<string> V = parsing::FASTASuffixes;
//		for (string s : V) {
//			cout << " " << s;
//		}
//		cout << endl;
//		cout << "NEXUS accepted suffixes: ";
//		V = parsing::NEXUSSuffixes;
//		for (string s : V) {
//			cout << " " << s;
//		}
//		cout << endl
//	);
	if (matchesIgnoreCase(suffix, parsing::FASTASuffixes)) {
		// parse as FASTA
		DEBUG(cout << "parsing " << filename << " as FASTA format file" << endl);
		parser = new parsing::FASTAParser(filename, this);
		ft = fastaFile;
	} else if (matchesIgnoreCase(suffix, parsing::NEXUSSuffixes)) {
		// parse as NEXUS
		DEBUG(cout << "parsing " << filename << " as NEXUSformat file" << endl);
		parser = new parsing::NEXUSParser(filename, this);
		ft = nexusFile;
	} else {
		stringstream ss;
		ss << "Cannot recognise the file name \'" << filename << "\': giving up.";
		throw new app_exception(ss.str());
	}
	DEBUG(cout << "setting project" << endl);
	parser->setProject(this);
	DEBUG(cout << "project set; now beginning parsing" << endl);
	parser->parse();
	DEBUG(cout << "checking sequence length... "; cout.flush());
	checkSequenceLength();
	DEBUG(cout << "ok\n");
	DEBUG(cout << "compressing alignment... "; cout.flush());
	compressAlignment();
	DEBUG(cout << "ok\n");
	checkSequenceLength(true);
}

/**
 * The only thing that matters in the pattern count conversion is the sign of each count;
 * in theory we multiply 1s and -1s and then the result has the right sign, but in practice
 * we can just count the negative entries and check the parity at the end.
 * Hence this matrix is renamed negH, and H = -2(negH) + 1 (entry-wise).
 */

Pattern<unsigned int>& encode(Pattern<unsigned int>& result, const Pattern<char>& pat) {
	result.resize(pat.size());
	char ch;
	unsigned int ntax = pat.size();
	for (unsigned int i = 0; i < ntax; ++i) {
		ch = pat[i];
		if (Project::encoding.find(ch) == Project::encoding.end()) {
			cout << "cannot find \'" << ch << "\' in encoding!" << endl;
			throw new app_exception("Whoops.");
		}
		unsigned int a = Project::encoding.at(ch);
		result[i] = a;
	}
	return result;
}

} /* namespace flatbush */
