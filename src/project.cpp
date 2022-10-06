/*
 * project.cpp
 *
 *  Created on: 25 Jul 2016
 *      Author: mac
 */

#include <algorithm>
#include <chrono>
#include <ctime>
#include <exception>
#include <iomanip>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <sys/stat.h>
#include <unordered_map>
#include <unordered_set>
#include <valarray>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Jacobi>
#include "project.h"
#include "../utility/appexception.h"
#include "../utility/debugging.h"
#include "../utility/phylogeny.h"
#include "../utility/parser.h"
#include "../utility/FASTAParser.h"
#include "../utility/NEXUSParser.h"
#include "../utility/pattern.h"
#include "../utility/tools.h"
//#include "../utility/CompactSequence.h"
#include "Climb.h"

using namespace std;

extern bool _debugging;
extern bool _save_climbs;
extern bool _silent;
extern bool _show_all_details;
//extern bool _show_taxon_names;
extern bool _verbose;
extern ofstream out;

/**
 * XXX Put in a data filter:
 * 	remove all sites with any ambiguity codes: can only handle upper case ACGT
 * 	convert all lower case acgt to ACGT respectively
 * 	remove all sites with gaps
 */
extern const unsigned int kNumStates;

// This is a 4x4 representation of the Hadamard matrix with a 1 if the sign is negative and 0 otherwise.
const int negH[4][4] = { { 0, 0, 0, 0 }, { 0, 1, 0, 1 }, { 0, 0, 1, 1 }, { 0, 1, 1, 0 } };

extern unsigned short NucsAsBits[256];
extern unsigned int numStatesFromNucCode[16];

namespace flatbush {

size_t card(const leafset& l) {
	return (static_cast<size_t>(count(l.begin(), l.end(), true)));
//	size_t c = 0;
//	for (unsigned int i = 0; i < l.size(); ++i) {
//		if (l[i]) {
//			++c;
//		}
//	}
//	return c;
}

string strmLeafset(const leafset& vec) {
	stringstream ss;
	ss << "{ ";
//	leafset* leaves = &vec;
//	static leafset flipped;
//	if (2*count(vec.begin(), vec.end(), true) > ntax) {
//		flipped = Project::flipLeafset(vec);
//		leaves = &flipped;
//	}
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
	DEBUG(bool _debugging(true));
	if (vec.size() == 0) {
		DEBUG(cout << "strmLeafsetBinary(const leafset& vec): vec is empty!" << endl);
		throw new app_exception("Attempting to pring an empty leafset: this shouldn't be called.");
	}
	if (vec[0] == true) {
		for (unsigned int i = 0; i < vec.size(); ++i) {
			if (vec[i] == true) {
				ss << '1';
			} else {
				ss << '0';
			}
//			ss << (vec[i]) ? '1' : '0';
		}
//			if (vec[i] == true) { // did use reverse order, with vec[vec.size()-1-i], but real order fits with Filo.
//				ss << '1'; //ss << (vec[i]) ? '1' : '0';
//			} else {
//				ss << '0';
//			}
	} else {
		for (unsigned int i = 0; i < vec.size(); ++i) {
			if (vec[i] == true) {
				ss << '0';
			} else {
				ss << '1';
			}
//			ss << (vec[i]) ? '0' : '1';
		}
	}
	DEBUG(cout << ss.str() << endl);
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

Project::~Project() {
	for (auto iter : walks) {
		delete iter.second;
	}
	if (D != nullptr) {
		for (unsigned int i = 0; i < ntax; ++i) {
			delete [] D[i];
		}
		delete [] D;
	}
	if (A != nullptr) {
		delete [] A;
	}
}

Project::Project() : parser(nullptr), ft(nullFile), dType(alignmentData),
		ntax(0), sf(sf_parsimony), seqLength(0), gapChar('-'),
		missingChar('?'), numStates(4), _hasPatternWeights(false), perts(2), sampleSize(1000),
		outfmt(fmt_binary_csv), _doAllSplits(false), _doLandscapeAnalysis(false),
		_suppress_labels(false),
		D(nullptr), A(nullptr),
		graphColMaxScoreRed(0.0), graphColMaxScoreGreen(1.0), graphColMaxScoreBlue(0.0), graphColMaxScoreAlpha(0.0),
		graphColMinScoreRed(0.0), graphColMinScoreGreen(0.0), graphColMinScoreBlue(0.5), graphColMinScoreAlpha(0.25)
 {
	initialise();
}

Project::Project(std::string filename) : parser(nullptr), perts(2), sampleSize(1000), outfmt(fmt_binary_csv),
		_doAllSplits(false), _doLandscapeAnalysis(false) {
	// load the file with this filename, tokenize it, then parse & create all the defined objects
	read(filename);
//	os = *cout;
}

void Project::addKnownSplit(std::string s) {
	// add a split that is to be checked against, for situations where
	// we know the true tree.
	splitsAsStrings.insert(s);
}

//void Project::addRawSequence(string& id, string& seq) {
//	cout << "addRawSequence(" << id << "," << seq << ")" << endl;
////	if (seqLength == 0) {
////		seqLength = seq.length();
////	} else if (seq.length() != seqLength) {
////		throw new app_exception("Cannot add sequences of different length to alignment");
////	}
//	rawAlignment.insert(std::pair<string,string>(id, seq));
//	ntax = rawAlignment.size();
//}

void Project::addSequence(string id, string seq) {
	rawAlignment[id] = seq;
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
	splitsAsLeafsets.push_back(split);
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
	splitsAsLeafsets.push_back(split);
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
	if (splitsAsLeafsets.size() == 0) {
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
	if (_show_all_details) {
		cout << "Known splits:\n";
		for (auto iter = splitsAsLeafsets.begin(); iter != splitsAsLeafsets.end(); ++iter)
			cout << "\t" << strmLeafset(*iter) << endl;
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
	for (leafset left : splitsAsLeafsets) {
//		DEBUG(cout << "left = " << strmLeafset(left) << '\t');
		right = left;
		right.flip();
//		DEBUG(cout << "right = " << strmLeafset(right) << endl);
		putHalfSplitMask(rows, left);
		putHalfSplitMask(columns, right);

		Eigen::MatrixXd F(rows.size(), columns.size());
//		Eigen::SparseMatrix<double> SM(rows)
		double error = 0.0;
		if (_show_all_details) {
			cout << endl << "Subflattening for split " << setw(ntax) << strmLeafset(left)
					<< ":" << strmLeafset(right) << ", " << describeSplit(left) << endl;
		}
		unsigned int r = 0;
		unsigned int c = 0;
//		cout << setw(static_cast<int>(log(A.getSequenceLength()))+15);
		if (_show_all_details) {
			cout << "\t" << nullPattern << "\t";
			for (auto col : columns) {
				cout << setw(10-ntax) << col << "\t";
			}
			cout << endl;
		}
		for (auto row : rows) {
			c = 0;
			if (_show_all_details) {
				cout << "\t" << row << "\t";
			}
			for (auto col : columns) {
				combine(pat, row, col);
				encode(idx, pat);
				F(r, c) = countSignedSum(idx);
				F(r, c) /= getSequenceLength();
				if (_show_all_details) {
					cout << setw(11) << F(r,c) << "\t";
				}
				++c;
			}
			if (_show_all_details) {
				cout << endl;
			}
			++r;
		}
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(F);
		Eigen::MatrixXd svalues = svd.singularValues().transpose();
		if (_show_all_details) {
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
		if (_show_all_details) {
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
			if (splitsAsStrings.find(leftString) != splitsAsStrings.end()
					|| splitsAsStrings.find(rightString) != splitsAsStrings.end()) {
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
		if (_show_all_details) {
			cout << error << ", " << strmLeafset(left) << ":" << strmLeafset(right) << ", " << describeSplit(left) << endl;
		}
//		}

	}

}

float Project::calcMeanBetweenClusterDistance(const leafset& left) const {
	float m(0.0);
	unsigned int size = count(left.begin(), left.end(), true);
	unsigned int denominator = size *(ntax-size);
	for (unsigned int i = 0; i < ntax; ++i) {
		if (left[i]) {
			// accumulate all the distances to everything NOT in this leafset
			for (unsigned int j = 0; j < ntax; ++j) {
				if (!left[j]) {
					m += D[i][j];
				}
			}
		} else {
			// accumulate all the distances to everything IN the leafset
			for (unsigned int j = 0; j < ntax; ++j) {
				if (left[j]) {
					m += D[i][j];
				}
			}
		}
	}
	m /= denominator;
	return m;
}

float Project::calcMeanWithinClusterDistance(const leafset& left) const {
	float m(0.0);
	unsigned int size = count(left.begin(), left.end(), true);
	unsigned int denominator = static_cast<int>(size * (size-1) / 2);
	for (unsigned int i = 0; i < ntax; ++i) {
		if (left[i]) {
			for (unsigned int j = i+1; j < ntax; ++j) {
				if (left[j]) {
					m += D[i][j];
				}
			}
		}
	}
	m /= denominator;
	return m;
}

void Project::calcDistanceMatrix(DistanceCorrection corr) {
	/**
	 * Given the alignment, calculate a distance matrix D, using the distance correction method corr.
	 */
	bool _debugging(true);
	D = new float*[ntax];
	for (unsigned int i = 0; i < ntax; ++i) {
		D[i] = new float[ntax];
	}
	if (A == nullptr) {
		setUpAlignmentPointers();
	}
	switch (corr) {
	case dc_p:
		for (unsigned int i = 0; i < ntax; ++i) {
			D[i][i] = 0.0;
			for (unsigned int j = i+1; j < ntax; ++j) {
				unsigned int h = calcHammingDistance(A[i], A[j]);
				D[i][j] = 1.0 * h / seqLength;
				D[j][i] = D[i][j];
			}
		}
		break;
	case dc_JC:
		for (unsigned int i = 0; i < ntax; ++i) {
			D[i][i] = 0.0;
			for (unsigned int j = i+1; j < ntax; ++j) {
				double p = calcPDistance(A[i], A[j]);
				DEBUG(cout << "p[" << i << "," << j << "] = " << p << endl);
				if (p >= 0.75) {
					cerr << "Fatal error: Jukes-Cantor Distance is not defined between taxa " << i << " and " << j;
					cerr << ", as p-distance is " << p << endl;
					throw new app_exception("Numerical error: saturated sequences.");
				}
				D[i][j] = -0.75 * log(1 - p / 0.75 );
				D[j][i] = D[i][j];
			}
		}
		break;
	case dc_HKY:
		throw new app_exception("Sorry, HKY distance is not yet implemented.");
		break;
	case dc_F81:
		throw new app_exception("Sorry, F81 distance is not yet implemented.");
		break;
	case dc_F84:
		throw new app_exception("Sorry, F84 distance is not yet implemented.");
		break;
	case dc_LogDet:
		throw new app_exception("Sorry, LogDet distance is not yet implemented.");
//		for (unsigned int i = 0; i < ntax; ++i) {
//			for (unsigned int j = i+1; j < ntax; ++j) {
//				/*
//				 * calculate F, a matrix of proportions of each A-A, A-C, ..., T-T pairing between the sequences
//				 */
//				D[i][j] = calcLogDetDistance(A[i], A[j]);
//			}
//		}
		break;
	default:
		break;
	}
	DEBUG(
		cout << "Distance Matrix:" << endl;
		for (unsigned int i = 0; i < ntax; ++i) {
			cout << i << '\t';
			for (unsigned int j = 0; j < ntax; ++j) {
				cout << D[i][j] << '\t';
			}
			cout << endl;
		}
	);

}

double Project::calcPDistance(const string* s, const string* t) {
	return (1.0 * calcHammingDistance(s, t) / s->length());
}

unsigned int Project::calcHammingDistance(const string* s, const string* t) const {
	bool _debugging(true);
	unsigned int h(0);
	if (s->length() != t->length()) {
		throw new app_exception("cannot calculate Hamming distance on strings of different lengths.");
	}
	for (unsigned int i = 0; i < s->length(); ++i) {
		h += (s->c_str()[i] != t->c_str()[i]) ? 1 : 0;
	}
	DEBUG(cout << "Hamming distance between strings\n\t" << s->c_str() << "\nand\n\t"
			<< t->c_str() << "\n\t= " << h << endl);
	return h;
}


float Project::calcLogDetDistance(const string* s, const string* t) const {
	float ld(0.0);
	return ld;
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
		compressAlignment();	// looks like this is being called again here.
		countPatterns(); // puts counts into patternCounts
		break;
	case patternWeightData:
		// should check the pattern weights are, e.g., all non-negative?
		try {
			for (auto iter = patternWeights.begin(); iter != patternWeights.end(); ++iter) {
				if (iter->first.length() != ntax) {
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

void Project::calcSeqFreqs() {
	/**
	 * For each sequence, calculate the number of instances of each character state;
	 * So for DNA sequences, the number of As, Cs, Gs and Ts in each sequence.
	 * E.g. if taxon1 has sequence AAGCCT the seqFreq would be [2,2,1,1].
	 *
	 * COMPLEXITY: O(ntax * seqlength * numstates)
	 */
	for (auto seq = rawAlignment.begin(); seq != rawAlignment.end(); ++seq) {
		seqFreqs[seq->first] = vector<int>(numStates, 0);	// create the count vector for this sequence name
		for (unsigned int i = 0; i < numStates; ++i) {
			seqFreqs[seq->first][i] = count(seq->second.begin(), seq->second.end(), symbols[i]);
				// use the built-in standard algorithm "count" to count number of matches of state i
		}
	}
}


void Project::captureAlignmentSettings(string str) {
	/**
	 * Scan the string (probably a file name) that should have important information about the alignments input.
	 * Patterns are for number of taxa, generating model, sequence length, and scale factor (branch length or height).
	 */
	bool _debugging(true);

	// trim any leading path:
	size_t pos = str.rfind('/', str.length());
	if (pos != string::npos) {
		str = str.substr(pos+1,str.length());
	}

	// regular expressions for info in the file name:
	regex patNTax("-n([0-9]+)");	// matches "-n21-" to set ntax=21
	regex patSeqLength("-l([0-9]+)");
	regex patSubstModel("([A-Z]+[0-9]*)");
	regex patScale("-[bh]([0-9]+\\.[0-9]+)");
	regex patSampleSize("-ss([0-9]+)");
	regex patNotes("^([a-zA-Z]+)-");
	smatch m;

	if (regex_search(str, m, patNotes)) {
		DEBUG(cout << "notes: " << m.str(1) << endl);
		inputNotes = m.str(1);
	}
	if (regex_search(str, m, patNTax)) {
		DEBUG(cout << "number of taxa = " << m.str(1) << endl);
		ntax = atoi(m.str(1).c_str());
	}
	if (regex_search(str, m, patSeqLength)) {
		DEBUG(cout << "sequence length = " << m.str(1) << endl);
		seqLength = atoi(m.str(1).c_str());
	}
	if (regex_search(str, m, patSubstModel)) {
		DEBUG(cout << "Substitution model = " << m.str(1) << endl);
		substitutionModel = m.str(1);
	}
	if (regex_search(str, m, patScale)) {
		DEBUG(cout << "scale = " << m.str(1) << endl);
		scaleFactor = atof(m.str(1).c_str());
	}
	if (regex_search(str, m, patSampleSize)) {
		DEBUG(cout << "sample size = " << m.str(1) << endl);
		sampleSize = atof(m.str(1).c_str());
	}
}

void Project::checkSequenceLength(bool nothrow) {
	bool _debugging(false);
	if (rawAlignment.size() == 0) {
		return;	// nothing to check as there are no sequences
	}
	auto iter = rawAlignment.begin();
	unsigned int sl = iter->second.size();
	DEBUG(cout << "first sequence = " << iter->second << endl);
	++iter;
	while (iter != rawAlignment.end()) {
		if (iter->second.size() != sl) {
			stringstream msg("two sequences of different length: sequence \'");
			msg << iter->first << " has length " << iter->second.size() << " but " << sl << " was expected.";
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

leafset Project::chooseRandomSplit() {
	leafset half(ntax, false);
	unsigned int card;
	do {
		for (unsigned int t = 0; t < ntax; ++t) { // don't include the last taxon, labelled (n-1).
			half[t] = (drand48() < 0.5);
//			if (drand48() < 0.5) {
//				half[t] = true;
//			} else {
//				half[t] = false;
//			}
		}
		shortenSplit(half);
		card = getCard(half);
	} while (card <= 1 || card >= ntax-1);
	return half;
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
//	/**
//	 * Now convert the strings into CompactSequence objects:
//	NO, REALLY DON'T
//	 */
//	compactAlignment.clear();
//	for (auto iter : rawAlignment) {
//		compactAlignment.insert(pair<string, string>(iter.first, CompactSequence(iter.second)));
//	}
}

/**
 * Search the Alignment A to find all the names of all the sequences;
 * then put them in a map (in order);
 * then assign each an index, which will be used in pattern counting (and possibly elsewhere).
 */
void Project::combine(Pattern<char>& result, const Pattern<char> x, const Pattern<char> y) {
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

unsigned int Project::countMatchingSplits(unordered_map<string, float> trueSplits, set<string> comparisonSplits) {
	/**
	 * Count and return how many splits -- represented as strings -- are in both "trueSplits" and "comparisonSplits"
	 */
	unsigned int numMatches(0);
	for (auto spl : trueSplits) {
		numMatches += comparisonSplits.count(spl.first);
	}
	return numMatches;
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
	checkSequenceLength(true);
	DEBUG(cout << "seqLength = " << seqLength << endl);
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
		const string& p = iter->first;
//		cout << p << ", ";
		int term = 0;
		for (unsigned int i = 0; i < patLength; ++i) {
//			if (pat[i] > 3) {
//				throw new app_exception("pat element is out of range");
//			}
//			if (p[i] != 'C' && p[i] != 'A' && p[i] != 'T' && p[i] != 'G') {
//				stringstream sstr;
//				sstr << "p = " << p << "; p[" << i << "] = " << p[i] << " is not in {A,C,G,T}" << endl;
//				throw new app_exception(sstr.str());
//			}
//			DEBUG(if (encoding.at(p[i]) < 0 || encoding.at(p[i]) > 3) {
//				throw new app_exception("encoding is out of range [0,3]");
//			})
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
	for (auto it : rawAlignment) {
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

typedef std::chrono::milliseconds millisecond;
typedef std::chrono::duration<float> fsec;

void Project::doLandscapeAnalysis() {
	bool _debugging(true);
	auto wallclock_start = std::chrono::system_clock::now();	// "actual" time for noting when the program is run
	auto steady_start = std::chrono::steady_clock::now();	// "steady" clock won't ever be adjusted down so ideal for accurate calculating duration of run
	std::time_t start_time = std::chrono::system_clock::to_time_t(wallclock_start);

	leafset right;
	vector<Pattern<char>> rows;
	vector<Pattern<char>> columns;
	Pattern<char> pat(ntax, '.');
	Pattern<unsigned int> idx;
	Pattern<char> nullPattern(ntax, ' ');
	unsigned int stepsToSuccess(0);	// number of steps in the climb to optimal and correct splits
	walks.clear();
	unsigned int walkLength;
	unsigned int hitRate(0);
	double score(0.0);
	if (leafLabels.size() == 0) {
		gatherTaxonLabels();
	}
	calcPatternWeights();
	if (sf == sf_cluster_distance) {
		calcDistanceMatrix(dc_JC);	// NOTE that I'm using Jukes-Cantor by default for now. JUST TESTING!!
	}
	try {
		ofstream os;
		if (_save_climbs) {
			os.open(name + "-climbs.csv");
			if (!os.good()) {
				throw new app_exception("Cannot open file for landscape output.  Sorreee!");
			}
			os << "StartSplit,StartScore,StopSplit,StopScore,SplitSize,PathLength\n";
		}
		if (ntax > 16) {
			_doAllSplits = false;
		}

		if (_doAllSplits) {
			sampleSize = (1 << (ntax-1));
			sampleSize -= ntax + 1;
			DEBUG(cout << " sampleSize = " << sampleSize << endl);
			gatherAllNontrivialHalfSplits(_verbose);
//			DEBUG(
//				for (auto split : splitsAsLeafsets) {
//					cout << strmLeafset(split) << endl;
//				}
//			);
			DEBUG(cout << "landscape analysis: number of splits to examine = " << splitsAsLeafsets.size() << endl);
			for (leafset left : splitsAsLeafsets) {
//				DEBUG(
//					for (const auto dix : splitsAsLeafsets) {
//						cout << "debugging: " << strmLeafset(dix) << endl;
//					}
//				);
				DEBUG(cout << "* Next leaf set: " << strmLeafset(left) << endl);
				score = getSplitError(left);
				leafset finish;
//				cout << "Start: " << strmLeafset(left) << ',' << score << endl;
				if (walks.find(left) == walks.end()) {	// that is, it's NOT already been encountered in a walk
					DEBUG(cout << "New starting point in base loop: " << strmLeafset(left) << endl);
	//				steepestDescentFromHere(left, walkLength);
					recursiveSteepestDescent(left, score, finish, walkLength, os);
					DEBUG(cout << "FINISHED New starting point in base loop" << endl);
				} else {
					finish = walks.at(left)->finish;
				}
				if (trueBranches.count(strmLeafsetBinary(finish)) > 0) {
					DEBUG(cout << "true split found! " << strmLeafset(finish) << " is in original tree." << endl);
					stepsToSuccess += walkLength;	// keeps track of how many steps it takes to get to any *correct* split
					hitRate++;
				}
//				cout << strmLeafset(left) << ',' << score << ',' << walkLength << endl;
			}

		} else {
			DEBUG(cout << "landscape analysis: number of splits to examine = " << sampleSize << endl);
//			unordered_set<leafset> startingSplits;
			for (unsigned int i = 0; i < sampleSize; ++i) {
				leafset left;
//				do {
					left = chooseRandomSplit();
//				} while (startingSplits.count(left) > 0);
				score = getSplitError(left);
				leafset finish;
				if (walks.find(left) == walks.end()) {
					recursiveSteepestDescent(left, score, finish, walkLength, os);
				} else {
					finish = walks[left]->finish;
				}
				if (trueBranches.count(strmLeafsetBinary(finish)) > 0) {
					stepsToSuccess += walkLength;	// keeps track of how many steps it takes to get to any *correct* split
					hitRate++;
				}
			}

		}
//		if (ntax > 16 || (sampleSize < 0)) {
//			if (sampleSize == 0) {
//				sampleSize = min(1000, 1 << (ntax-1) - ntax + 1);
//			}
//			_doAllSplits = false;
//			// do it with random starts,  just 1000 say
//		} else {
//
//		}
		if (_save_climbs) {
			os.close();
		}

		if (_verbose) {
			cout << "Domains of attraction:\n";
		}
		os.open(name + "-doa.csv");
		if (!os.good()) {
			throw new app_exception("Cannot open file for domains of attraction. Drat.");
		}
		os << "Split,SplitSize,FinalScore,DoASize" << endl;
		if (_verbose) {
			cout << "Split,SplitSize,FinalScore,DoASize" << endl;
		}
		unsigned int sumDoASize(0);
		multimap<int, leafset> doas;
		for (auto pr : domainOfAttraction) {
			doas.insert(pair<int, leafset>(pr.second, pr.first));
		}
//		doas.insert(domainOfAttraction.begin(), domainOfAttraction.end());
			// inserting the lot from the unordered map into the order map, so the output is nice.
		string bsplit;
		set<string> LOpt;	// the set of all locally optimal splits, as binary strings
		for (auto iter = doas.rbegin(); iter != doas.rend(); ++iter) {
			bsplit = strmLeafsetBinary(iter->second);
			os << bsplit
					<< ',' << std::count(iter->second.begin(), iter->second.end(), 1)
					<< ',' << getSplitError(iter->second)
					<< ',' << iter->first
					<< endl;
			if (_verbose) {
				cout << bsplit
							<< ',' << std::count(iter->second.begin(), iter->second.end(), 1)
							<< ',' << getSplitError(iter->second)
							<< ',' << iter->first
							<< endl;
			}
			LOpt.insert(bsplit);
			sumDoASize += iter->first;
		}
		os.close();

		bool _save_summary_statistics(true);
		string flatbushSummaryFile("flatbush-summary.csv");
		struct stat buf;
		if (stat(flatbushSummaryFile.c_str(), &buf) == -1) {
			os.open(flatbushSummaryFile);
			os << "date,note,model,n,l,scale,splitScoreMethod,sampleSize,nLOptSplits,nCorrectLOptSplits,hitRate,meanStepsToSuccess,elapsed\n";
		} else {
			os.open(flatbushSummaryFile, std::ofstream::app);	// open a standard file for appending these results
		}

		auto steady_end = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration<float>(steady_end - steady_start);
		millisecond duration = std::chrono::duration_cast<millisecond>(elapsed);

		if (_save_summary_statistics) {
//			ostringstream oss;
//			oss <<
			struct std::tm* ptm = std::localtime(&start_time);
			os << put_time(ptm, "\"%Y-%m-%d-%H:%M:%S\"") << ',';
			os << inputNotes << ',';
			os << substitutionModel << ',';
			os << ntax << ',';
			os << seqLength << ',';
			os << scaleFactor << ',';
			os << strmSplitScoreMethod() << ',';
			os << sampleSize << ',';
			os << doas.size() << ',';
			os << countMatchingSplits(trueBranches, LOpt) << ',';
			os << 1.0 * hitRate / sampleSize << ',';
			os << 1.0 * stepsToSuccess / hitRate << ',';
			os << elapsed.count() << endl;
			os.close();
		}


		if (_verbose) {
			cout << "Sum of domain sizes = " << sumDoASize << endl;
			cout << "Landscape analysis complete." << endl;
		}
		exportLandscape();
	} catch (app_exception& e) {
		cerr << "Exception in Project::doLandscapeAnalysis()" << endl;
	}
}

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

void Project::exportLandscape() {
	/**
	 * Create a dot file (for GraphViz) with each node being a split,
	 * 	edges between splits if they're adjacent under the current perturbation,
	 * 	and a weight for each split corresponding to its score.
	 * Not sure how to convey weights in dot... except by colour maybe?
	 */
	ofstream out;
	try {
		string landscapeFilename = name + "-graph.dot";
		out.open(landscapeFilename);
		if (!out.good()) {
			throw new app_exception("Cannot open dot file for exporting the landscape graph.");
		}
		double maxError = -1e20;	// XXX should make this a max double whatever it is
		double minError = 1e20;
		out << "graph " << "Landscape" << " {\n"
				<< "\tlayout=\"neato\"\n"
				<< "\tnode [shape=\"none\", height=0 style=\"filled\" fontcolor=\"black\"];\n"
				<< "\tedge [penwidth=2, color=\"lightgrey\"];\n"
				<< "\tgraph [nodesep=0.5, pad=\"0.5\"];\n";
		if (_suppress_labels) {
			out << "\tnode [shape=\"circle\", height=0.5]";	// a bit clunky just overriding the previous command but hey.
		}
		for (auto splt : splitError) {
			const leafset& a = splt.first;
			double err = getSplitError(a);
			maxError = max(maxError, err);
			minError = min(minError, err);
		}
//		double minHue(0.2);
//		double hueRange(0.2);
		out << "\t#[maxError = " << maxError << ']' << endl;
		out << "\t#[minError = " << minError << ']' << endl;

		out << "\tKEY [shape=\"none\" fontsize=18 label=<\n";	// XXX magic font size
		out << "\t\t<table cellpadding=\"1\" cellborder=\"0\" cellspacing=\"0\" border=\"0\">\n";
		out << "\t\t<tr align=\"center\" ><td><b>KEY</b></td></tr>" << endl;
		for (int i = 0; i < 11; ++i) {
			double val = 0.1 * i;
			string colstr("#");
			colstr += unitTo256Hex(graphColMinScoreRed + (graphColMaxScoreRed - graphColMinScoreRed)*val);	// red
			colstr += unitTo256Hex(graphColMinScoreGreen + (graphColMaxScoreGreen - graphColMinScoreGreen)*val);	// green
			colstr += unitTo256Hex(graphColMinScoreBlue + (graphColMaxScoreBlue - graphColMinScoreBlue)*val);	// blue
			colstr += unitTo256Hex(graphColMinScoreAlpha + (graphColMaxScoreAlpha - graphColMinScoreAlpha)*val);	// alpha
			if (i == 0) {
				out << "\t\t<tr><td bgcolor=\"" << colstr << "\">min(best)</td></tr>\n";
			} else if (i == 10) {
				out << "\t\t<tr><td bgcolor=\"" << colstr << "\">max(worst)</td></tr>\n";
			} else {
				out << "\t\t<tr><td bgcolor=\"" << colstr << "\">min+" << i << "0%</td></tr>\n";
			}
		}
		out << "\t\t</table>\n\t>]\n";

//		// KEY subgraph:
//		out << "\tsubgraph cluster0 {" << endl;
//		out << "\t\tlabel = \"Key\";" << endl;
//		out << "\t\tstyle = filled;" << endl;
//		out << "\t\tcolor = cornsilk;" << endl;
//		out << "\t\tnode [shape=\"none\"];" << endl;
//		// write out 11 nodes at 10% intervals from minScore to maxScore as a key
//		for (unsigned int i = 0; i < 11; ++i) {
//			if (i == 0) {
//				out << "\t\tk" << i << "0" << " [label=\"min(best)" << "\"";
//			} else if (i == 10) {
//				out << "\t\tk" << i << "0" << " [label=\"max(worst)" << "\"";
//			} else {
//				out << "\t\tk" << i << "0" << " [label=\"min+" << i << "0%" << "\"";
//			}
//			out << ", style=\"filled\"";
//			out << ", fillcolor=\"#";
//			double val = 0.1*i;
//			// output red green blue alpha
//			out << unitTo256Hex(graphColMinScoreRed + (graphColMaxScoreRed - graphColMinScoreRed)*val);	// red
//			out << unitTo256Hex(graphColMinScoreGreen + (graphColMaxScoreGreen - graphColMinScoreGreen)*val);	// green
//			out << unitTo256Hex(graphColMinScoreBlue + (graphColMaxScoreBlue - graphColMinScoreBlue)*val);	// blue
//			out << unitTo256Hex(graphColMinScoreAlpha + (graphColMaxScoreAlpha - graphColMinScoreAlpha)*val);	// alpha
//
//			// leftover from using HSV, which I might go back to including:
////			out << ' ';
////			out << (1.0-val)*(1.0-val);	// how intense (saturated) is the colour: very pale through to bright and saturated.
////			out << ' ';
////			out << 0.8;	// this
//
////			out << val;
////			out << minHue + val;
////			if (graphColHue > 0.0) {
////				out << minHue + hueRange*pow(val/hueRange, graphColHue);
////			} else {
////				out << minHue + hueRange*pow((1-val/hueRange), -graphColHue);
////			}
////			out << 0.5;
////			if (graphColSaturation > 0.0) {
////				out << pow(val, graphColSaturation);
////			} else {
////				out << pow((1-val), -graphColSaturation);
////			}
////			if (graphColValue > 0.0) {
////				out << pow(val, graphColValue);
////			} else {
////				out << pow((1-val), -graphColValue);
////			}
//
//			out << "\", fontcolor=\"black";
////			out << "\", fontcolor=\"0.2 " << (1-val) << ' ' << (1-val);
//			out << "\"];" << endl;
//		}
//		for (unsigned int i = 0; i < 10; ++i) {
//			out << "\t\tk" << i << "0 -- k" << (i+1) << "0 [len=0];\n";
//		}
//		out << "\t\tk100 -- k00 [style=\"invis\"];\n";
//		out << "\t}" << endl;
//		// end KEY

//		string bsplit;
//		for (auto spl: splitError) {
//			const leafset& a = spl.first;
////			cout << "Getting split error for " << strmLeafsetBinary(a) << endl;
//			double val = hueRange*(getSplitError(a)-minError)/(maxError-minError);
//			stringstream label;
//			label << strmLeafset(a);
//			if (domainOfAttraction[a] > 0) {
//				label << " [" << domainOfAttraction[a] << ']';
//			}
//			out << "\ts" << strmLeafsetBinary(a) << " [label=\"" << label.str() << "\"";
//			out << ", style=\"filled\"";
//			if (trueBranches.count(strmLeafsetBinary(a))) {
//				out << ", shape=\"square\"";
//				out << ", color=\"orange";
//			} else {
//				out << ", color=\"";
//				out << " 0.4 " << (1.0-val)*(1.0-val) << " 0.8 ";
////				if (graphColHue > 0.0) {
////					out << minHue + hueRange*pow(val/hueRange, graphColHue);
////				} else {
////					out << minHue + hueRange*pow((1-val/hueRange), -graphColHue);
////				}
////				if (graphColSaturation > 0.0) {
////					out << pow(val, graphColSaturation);
////				} else {
////					out << pow((1-val), -graphColSaturation);
////				}
////				if (graphColValue > 0.0) {
////					out << pow(val, graphColValue);
////				} else {
////					out << pow((1-val), -graphColValue);
////				}
//			}
//			out << "\", fontcolor=\"black\"];";
//			out << " // score=" << getSplitError(a);
//			out << endl;
//
////			out << "\", fontcolor=\"0.2 " << (1-val) << ' ' << (1-val);
////			out << "\"];" << endl;
//		}
		bool _showOnlyClimbEdges(true);
		if (_showOnlyClimbEdges) {
			for (auto e : climbEdges) {

				const leafset a(e.first);
//				double val = hueRange*(getSplitError(a)-minError)/(maxError-minError);
				double val = (getSplitError(a)-minError)/(maxError-minError);
				stringstream label;
				label << strmLeafset(a);
				out << "\ts" << strmLeafsetBinary(a) << " [label=\"";
				if (_suppress_labels) {
					out << " ";
				} else {
					out << label.str();
				}
				out << "\"";
//				out << ", style=\"filled\"";
				if (trueBranches.count(strmLeafsetBinary(a))) {
					out << ", shape=\"rect\"";
					out << ", color=\"orange";
				} else {
					out << ", color=\"#";
//					out << " 0.4 " << (1.0-val)*(1.0-val) << " 0.8 ";
					out << unitTo256Hex(graphColMinScoreRed + (graphColMaxScoreRed - graphColMinScoreRed)*val);	// red
					out << unitTo256Hex(graphColMinScoreGreen + (graphColMaxScoreGreen - graphColMinScoreGreen)*val);	// green
					out << unitTo256Hex(graphColMinScoreBlue + (graphColMaxScoreBlue - graphColMinScoreBlue)*val);	// blue
					out << unitTo256Hex(graphColMinScoreAlpha + (graphColMaxScoreAlpha - graphColMinScoreAlpha)*val);	// alpha
				}
//				out << "\", fontcolor=\"black\"];";
				out << "\"];";
				out << " // score=" << getSplitError(a);
				out << endl;

				const leafset b = e.second;
				if (domainOfAttraction[b] > 0) {
//					val = hueRange*(getSplitError(b)-minError)/(maxError-minError);
					val = (getSplitError(b)-minError)/(maxError-minError);
					stringstream label;
					label << strmLeafset(b);
					out << "\ts" << strmLeafsetBinary(b) << " [label=\"";
					if (_suppress_labels) {
						out << " ";
					} else {
						out << label.str();
					}
					out << "\"";
//					out << ", style=\"filled\"";
					if (trueBranches.count(strmLeafsetBinary(b))) {
						out << ", shape=\"square\"";
						out << ", color=\"orange";
					} else {
						out << ", color=\"#";
//						out << " 0.4 " << (1.0-val)*(1.0-val) << " 0.8 ";
						out << unitTo256Hex(graphColMinScoreRed + (graphColMaxScoreRed - graphColMinScoreRed)*val);	// red
						out << unitTo256Hex(graphColMinScoreGreen + (graphColMaxScoreGreen - graphColMinScoreGreen)*val);	// green
						out << unitTo256Hex(graphColMinScoreBlue + (graphColMaxScoreBlue - graphColMinScoreBlue)*val);	// blue
						out << unitTo256Hex(graphColMinScoreAlpha + (graphColMaxScoreAlpha - graphColMinScoreAlpha)*val);	// alpha
					}
//					out << "\", fontcolor=\"black\"";
					out << "\"];";
					out << " // score=" << getSplitError(b);
					out << endl;
				}
				if (e.first != e.second) {
					out << "\ts" << strmLeafsetBinary(e.first) << " -- s" << strmLeafsetBinary(e.second) << ";\n";
				}
			}
		} else {
			for (vector<leafset>::const_iterator i = splitsAsLeafsets.begin(); i != splitsAsLeafsets.end(); ++i) {
				for (vector<leafset>::const_iterator j = i+1; j!= splitsAsLeafsets.end(); ++j) {
					if (HammingDistance(*i, *j) == 2) {
						out << "\ts" << strmLeafsetBinary(*i) << " -- s" << strmLeafsetBinary(*j) << ";" << endl;
					}
				}
			}
		}
		out << "}\n";
		out.close();

		string splitsFilename = name + "-splits.csv";
		out.open(splitsFilename);
		if (!out.good()) {
			throw new app_exception("Cannot open dot file for exporting the landscape graph.");
		}
//		for (auto spl : trueBranches) {
//			cout << "---" << spl.first << ":" << spl.second << endl;
//		}
		out << "bsplit,leafset,size,score" << endl;
		for (const leafset& a : splitsAsLeafsets) {
			stringstream label;
			string bsplit(strmLeafsetBinary(a));
//			cout << bsplit << endl;
			label << strmLeafset(a);
			out << bsplit << ", \"" << label.str() << "\"," << getCard(a) << ","
					<< getSplitError(a) << "," << trueBranches.count(bsplit) << endl;
		}
		out.close();

	}
	catch (app_exception& e) {
		cerr << "exporting landscape failed." << e.what() << endl;
		return;
	}
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
	climbEdges.push_back(make_pair(L, bestNeighbour));
	return bestScore;
}

leafset Project::flipLeafset(const leafset& L) {
/**
 * Make and return a leafset that is the "flipped" description of the original.
 */
	leafset F(L);
	for (unsigned int i = 0; i < L.size(); ++i) {
		F[i] = !F[i];
	}
	return F;
}

bool Project::hasSequence(string seqID) {
	return (rawAlignment.find(seqID) != rawAlignment.end());
}

void Project::initialise() {
	char buffer[20];
	struct tm *timeinfo;
	time_t rawtime;
	time (&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(buffer, 20, "%Y-%m-%d-%H-%M-%S", timeinfo);
	string timeStr = buffer;
	name = "flatbush-project-" + timeStr;
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

// Alphabetical order of character states A C G T by default:
//but just out of interest....
void Project::gatherAllNontrivialHalfSplits(bool _verbose) {
	bool _debugging(false);
//	int ntax = A.getNumTaxa();
	if (ntax > 20) {
		throw new app_exception("that's too many splits, sorry. Stick to 20 taxa or fewer for this method.");
	}
	splitsAsLeafsets.clear();	// empties the vector and sets its size as 0
	leafset half(ntax, false);
	for (int i = 3; i < (1<<(ntax-1)); ++i) { // XXX At some point a Gray code would be cool to implement for this.
		// check out coding answer here: https://stackoverflow.com/questions/17490431/gray-code-increment-function
		unsigned int card = 0;
		for (unsigned int t = 0; t < ntax; ++t) { // don't include the last taxon, labelled (n-1).
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
			for (unsigned int x = 0; x < ntax; ++x) {
				if (half[x]) {
					cout << x << " ";
				}
			}
			cout << "}" << endl
		);
		shortenSplit(half);
//		if (static_cast<unsigned int>(2*count(half.begin(), half.end(), true)) > ntax) {
//			half = flipLeafset(half);
//		}
		splitsAsLeafsets.push_back(half);
//		if (_verbose) {
//			cout << strmLeafset(half) << " ";
//		}
	}
	cout.flush();
}

void Project::gatherTaxonLabels() {
	for (auto taxSeq : rawAlignment) {
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

unsigned int Project::getCard(const leafset& leaves) {
	unsigned int card(0);
	for (auto b = leaves.begin(); b != leaves.end(); ++b) {
		card += (*b) ? 1 : 0;
	}
	return card;
}

leafset Project::getRandomSplit() {
	leafset L(ntax);
	do {
		for (unsigned int i = 0; i < ntax; ++i) {
			L[i] = (drand48() < 0.5) ? 1 : 0;
		}
		shortenSplit(L);
	} while (isTrivialSplit(L));
	return L;
}

string& Project::getSequence(unsigned int seqID) {
	return *(A[seqID]);
}
const string& Project::getSequence(unsigned int seqID) const {
	return *(A[seqID]);
}
string& Project::getSequence(const string& seqID) {
	return rawAlignment.at(seqID);
}
const string& Project::getSequence(const string& seqID) const {
	return rawAlignment.at(seqID);
}

double Project::getSplitError(const leafset& split) {
	bool _debugging(false);
	DEBUG(cout << "Project::getSplitError(" << strmLeafset(split) << "}:" << endl);
	switch (sf) {
	case sf_cluster_distance:
		return getSplitErrorClustering(split);
	case sf_singular_value:
		return getSplitErrorSVD(split);
	case sf_parsimony:
		return getSplitErrorParsimony(split);
	case sf_state_distribution:
		return getSplitErrorStateDistribution(split);
	default:
		return 1.0;
	}
}

double Project::getSplitErrorClustering(const leafset& left) {
	/**
	 * This split score is the difference in mean distance of taxa across the split vs within each part.
	 */
	bool _debugging(true);
	if (D == nullptr) {
		throw new app_exception("Cannot perform the clustering split score with no distance matrix defined.");
	}
	if (splitError.find(left) != splitError.end()) {
		return splitError.at(left);
	}
	DEBUG(cout << "getSplitErrorClustering method" << endl);
	double P(0.0);
	leafset right(left);
	right.flip();
	float dA, dB, dAB;
	dA = calcMeanWithinClusterDistance(left);
	dB = calcMeanWithinClusterDistance(right);
	dAB = calcMeanBetweenClusterDistance(left);
	DEBUG(cout << "dA = " << dA << "; dB = " << dB << "; dAB = " << dAB << endl);
	P = 0.5*(dA + dB) - dAB;
	splitError[left] = P;
	return P;
}


double Project::getSplitErrorStateDistribution(const leafset& left) {
	/**
	 * This split score is the distance between the two vectors of base frequencies
	 * on each side of the split.
	 */
	static bool _seqFreqDistributionKnown(false);

	if (!_seqFreqDistributionKnown) {
		calcSeqFreqs();
		_seqFreqDistributionKnown = true;
	}
	bool _debugging(false);
	if (left.size() == 0) {
		throw new app_exception("leafset is empty in Project::getSplitErrorSVD()");
	}
	if (splitError.find(left) != splitError.end()) {
		return splitError.at(left);
	}
	bool _newway(true);
	double P(0);
	if (_newway) {
		valarray<float> leftCharacterProportion(numStates), rightCharacterProportion(numStates);
		leftCharacterProportion.resize(numStates);	// [0.0,0.0,0.0,0.0] for DNA (currently hard-wired)
		// sum the symbol/base frequency counts for the complete sequence on each side of the split; don't do the individual sites.
		for (auto taxon = leafLabels.begin(); taxon != leafLabels.end(); ++taxon) { // taxon is (string, int) and maps taxon names to their index
			if (left[taxon->second]) {
				for (unsigned int i = 0; i < numStates; ++i) {
					leftCharacterProportion[i] += seqFreqs[taxon->first][i];
				}
			} else {
				for (unsigned int i = 0; i < numStates; ++i) {
					rightCharacterProportion[i] -= seqFreqs[taxon->first][i];
				}
			}
		}
		unsigned int leftSize(std::count(left.begin(), left.end(), 1));
		leftCharacterProportion /= (leftSize * seqLength);
		unsigned int rightSize(ntax - leftSize);
		rightCharacterProportion /= (rightSize *seqLength);
		float diff(1.0 * numStates);
		for (unsigned int i = 0; i < numStates; ++i) {
			diff -= (leftCharacterProportion[i] - rightCharacterProportion[i]) * (leftCharacterProportion[i] - rightCharacterProportion[i]);
		}
		P = sqrt(diff);
	} else {
		map<char, int> diffV;
		int patternScore(0);
		int diff(0);
		for (auto pat = patternCounts.begin(); pat != patternCounts.end(); ++pat) {
			// accumulate difference between counts of states on left and right of split:
			for (char state : getStates()) {
				diffV[state] = 0;
			}
			for (unsigned int taxon = 0; taxon < ntax; ++taxon) {	// for each taxon
				if (left[taxon]) {											// if taxon is on left side of split
					diffV[pat->first[taxon]] += 1;						// increase count of this state
				} else {
					diffV[pat->first[taxon]] -= 1;						// else decrement
				}
			}
			// sum (leftV[i]-rightV[i])^2
			patternScore = 0;
			unsigned int size(std::count(left.begin(), left.end(), 1));
			diff = size*size + (ntax-size)*(ntax-size);	// |left|^2 + |right|^2;
				// this is the maximum possible score for a split of this size
			for (char state : getStates()) {
				DEBUG(cout << "diffV[" << state << "] = " << diffV[state] << endl);
				diff -= diffV[state]*diffV[state];
			}
			patternScore += diff * pat->second;
			DEBUG(cout << "contribution = sqrt(sum diff^2)*count = sqrt("
					<< diff << ")*" << pat->second << " = " << (sqrt(diff) *pat->second)
					<< endl);
			DEBUG(cout << "State distribution split score = " << patternScore << endl);
			P += patternScore;
		}
		P /= static_cast<double>(seqLength);
	}
	splitError[left] = P;
	return P;
}

double Project::getSplitErrorSVD(const leafset& left) {
	bool _debugging(false);
	if (left.size() == 0) {
		throw new app_exception("leafset is empty in Project::getSplitErrorSVD()");
	}
	if (splitError.find(left) != splitError.end()) {
		return splitError.at(left);
	}
//	bool _verbose(true);
	DEBUG(cout << "** input split (left) = " << strmLeafset(left) << '\t');
	leafset right(left);
	right.flip();
	DEBUG(cout << "right = " << strmLeafset(right) << endl);
	vector<Pattern<char>> rows;
	vector<Pattern<char>> columns;
	Pattern<char> pat(ntax, '.');
	Pattern<unsigned int> idx;
	putHalfSplitMask(rows, left);
	putHalfSplitMask(columns, right);

	Eigen::MatrixXd F(rows.size(), columns.size());
//		Eigen::SparseMatrix<double> SM(rows)
	double error = 0.0;
	if (_show_all_details) {
		cout << endl << "Subflattening for split " << setw(ntax) << strmLeafset(left)
				<< ":" << strmLeafset(right) << ", " << describeSplit(left) << endl;
	}
	unsigned int r = 0;
	unsigned int c = 0;
//		cout << setw(static_cast<int>(log(A.getSequenceLength()))+15);
	if (_show_all_details) {
		Pattern<char> nullPattern(ntax, ' ');
		cout << "\t" << nullPattern << "\t";
		for (auto col : columns) {
			cout << setw(10-ntax) << col << "\t";
		}
		cout << endl;
	}
	for (auto row : rows) {
		c = 0;
		if (_show_all_details) {
			cout << "\t" << row << "\t";
		}
		for (auto col : columns) {
			combine(pat, row, col);
			encode(idx, pat);
			F(r, c) = countSignedSum(idx);
			F(r, c) /= getSequenceLength();
			if (_show_all_details) {
				cout << setw(11) << F(r,c) << "\t";
			}
			++c;
		}
		if (_show_all_details) {
			cout << endl;
		}
		++r;
	}
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(F);
	Eigen::MatrixXd svalues = svd.singularValues().transpose();
	if (_show_all_details) {
		cout << "\tSingular Values found by Jacobi SVD: " << svalues << endl;
	}
	/**
	 * ERROR calculation:
	 * Sum from (k) to (the number of columns of F) the square of the singular value
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
	if (_show_all_details) {
		cout << "\tError from rank " << numStates << " : " << error << endl;
	}
	if (out.good()) {
//			out.precision(6);
//			out.setf(std::ios::fixed, std:: ios::floatfield);
		// TODO In this part, output a flag too, to say whether the split is in the tree.
//			DEBUG(cout << idx << "\t" << error << endl);
		string leftString = strmLeafsetBinary(left);
		string rightString = strmLeafsetBinary(right);
//		DEBUG(cout << setw(8) << error << ", " << leftString << ", " << rightString << endl);
//			cout << "split: " << leftString << "\t" << rightString << endl;
		int _inTree = 0;
		if (splitsAsStrings.find(leftString) != splitsAsStrings.end()
				|| splitsAsStrings.find(rightString) != splitsAsStrings.end()) {
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

double Project::getSplitErrorParsimony(const leafset& left) {
	unsigned int P(0);
	/**
	 * Algorithm:
	 * For each site i (i=1,...,numPatterns)
	 * 	Let L be the set of characters at site i at the leaves given by leafset left
	 * 	Let R be the set of characters at site i at the leaves given by leafset ([n] - left)
	 * 	P += |L| + |R| + (1 if L \intersect R is empty)
	 * 	XXX I'm not convinced that this algorithm is correct!
	 */
	if (splitError.find(left) != splitError.end()) {
		return splitError.at(left);
	}
//	DEBUG(cout << "input split (left) = " << strmLeafset(left) << '\t');
	bool _debugging(false);
//	leafset right(left);
//	right.flip();
	unsigned int L(0);
	unsigned int R(0);
	unsigned int patternScore(0);
	for (auto pat = patternCounts.begin(); pat != patternCounts.end(); ++pat) {
		L = 0;
		R = 0;
		for (unsigned int i = 0; i < ntax; ++i) {
			if (left[i]) {
				L |= NucsAsBits[static_cast<unsigned int>(pat->first[i])];
			} else {
				R |= NucsAsBits[static_cast<unsigned int>(pat->first[i])];
			}
		}
		patternScore = numStatesFromNucCode[L] + numStatesFromNucCode[R] - 2;
		patternScore += (L & R) ? 0 : 1;	// add 1 if L and R have no states in common
//		if (pat->second != 1) {
//			cout << "This one has more than one instance!\n";
//		}
		patternScore *= pat->second;
		DEBUG(cout << "Parsimony pattern score = " << patternScore << endl);
		P += patternScore;
	}
	splitError[left] = P / static_cast<double>(seqLength);
	return P;
}

vector<leafset>* Project::getSplits() {
	return &splitsAsLeafsets;
}

unsigned int Project::HammingDistance(const leafset& a, const leafset& b) {
	/**
	 * Count the number of differences between a and b.
	 */
	bool _debugging(false);
	unsigned int h(0);
	for (unsigned int i = 0; i < a.size(); ++i) {
		h += (a[i] != b[i]) ? 1 : 0 ;
	}
	DEBUG(cout << "HammingDistance(" << strmLeafset(a) << ',' << strmLeafset(b) << ") = " << h << endl);
	return h;
}

bool Project::isTrivialSplit(const leafset& split) {
	unsigned int size = card(split);
	bool _debugging(false);
	bool _trivial((size < 2) || (size > ntax-2));
	DEBUG(cout << "testing for trivial split: split = " << strmLeafset(split)
			<< "; size = " << size << "; trivial = " << _trivial << endl);
	return _trivial;
}

void Project::putAllAdjacentSplits(vector<leafset>& N, const leafset L) {
	bool _debugging(false);
	N.clear();
	if (perts & 1) { // perturbation #1
		// First perturbation: move a single taxon from one side to the other:
		for (unsigned int i = 0; i < ntax; ++i) {
			leafset a(L);
			a[i] = 1 - a[i]; // flip the ith bit
			if (isTrivialSplit(a)) {
				continue;
			}
			DEBUG(cout << "original split: " << strmLeafset(L) << "; next split by single flip = " << strmLeafset(a) << endl);
			shortenSplit(a);
			N.push_back(a);
		}
	}
	if (perts & 2) {
		// Second perturbation: swap two taxa:
		for (unsigned int i = 0; i < ntax; ++i) {
	//		if (L[i] == 0) {
	//			continue;
	//		}
			for (unsigned int j = i+1; j < ntax; ++j) {
				if (L[i] == L[j])
					continue;
				leafset a(L);
				a[i] = 1-a[i];
				a[j] = 1-a[j];
				DEBUG(cout << "original split: " << strmLeafset(L) << "; next split by double flip = " << strmLeafset(a) << endl);
	//			if (isTrivialSplit(a)) {
	//				continue;
	//			}
	//			if (a[ntax-1]) {
	//				// replace with the complement as this set includes (ntax-1)
	//				for (unsigned int pos = 0; pos < a.size(); ++pos) {
	//					a[pos] = 1 - a[pos];
	//				}
	//			}
				shortenSplit(a);
				N.push_back(a);
			}
		}
	}
	DEBUG(
		cout << "All neighbouring splits of " << strmLeafset(L) << ":\n";
		for (leafset& s : N) {
			cout << '\t' << strmLeafset(s) << endl;
		}
	)
}

double Project::recursiveSteepestDescent(const leafset& currentSplit, double currentScore, leafset& finish,
		unsigned int& walkLength, ostream& os) {
	/**
	 * Check all neighbours for that split with the minimal score
	 * Select that split, b
	 * If walk(b) does not exist then
	 * 	recursiveSteepestDescent(b)
	 * walk(current) to contain walks(b).finish(), walks(b).finalScore(), and walks(b).walkLength+1
	 */
	bool _debugging(true);
	leafset b;
	DEBUG(cout << "recursiveSteepestDescent: " << strmLeafset(currentSplit) << endl);
//	splitsAsLeafsets.push_back(currentSplit); // TODO XXX THIS IS THE BIT THAT'S FUCKING THINGS UP
	double bestScore = findBestNeighbouringSplit(b, currentSplit, getSplitError(currentSplit));
	DEBUG(cout << "bestScore = " << bestScore << endl);
	if (walks.find(b) != walks.end()) {
		// we have this descent / climb already from b to a local opt, so just tack this on to the currentSplit climb:
		Climb* desc = walks.at(b);
		DEBUG(cout << "This is a known walk from " << strmLeafset(b) << " to " << strmLeafset(desc->finish)
				<< "; final score = " << desc->finalScore << " and length " << desc->walkLength << ".\n");
		Climb* climb = new Climb(currentSplit, getSplitError(currentSplit), desc->finish, desc->walkLength+1, desc->finalScore);
		DEBUG(cout << "New climb: " << strmLeafset(currentSplit) << " -> " << strmLeafset(climb->finish)
				<< ", score = " << climb->finalScore << " in " << climb->walkLength << " steps.\n");
		walks[currentSplit] = climb;
		walkLength = desc->walkLength+1;
		finish = desc->finish;
		if (os.good()) {
			os << (*climb) << endl;
		}
		if (_verbose) {
			cout << (*climb) << endl;
		}
		DEBUG(cout << "size of walks map = " << walks.size() << endl);
		if (finish.size() == 0) {
			throw new app_exception("[A] \'finish\' leafset is empty");
		}
		domainOfAttraction[finish]++;
		DEBUG(cout << "Domain of attraction for split " << strmLeafset(finish) << " is increased to " << domainOfAttraction[finish] << endl); // XXXHERE
		return desc->finalScore;
	}
	// Else, we don't have this walk recorded already.
	DEBUG(cout << "Have not encountered this walk before.\n");
	DEBUG(cout << "bestScore = " << bestScore << "; currentScore = " << currentScore
			<< "; difference = " << (bestScore - currentScore) << endl);
	if (b == currentSplit) {
		// If there are no more improving neighbours, this is a walk of length 0
		Climb* climb = new Climb(currentSplit, getSplitError(currentSplit), currentSplit, 0, bestScore);
		DEBUG(cout << "there are no more improving neighbours; this is a walk of length 0" << endl);
		finish = currentSplit;
		walkLength = 0;
		walks[currentSplit] = climb;
		if (os.good()) {
			os << (*climb) << endl;
		}
		if (_verbose) {
			cout << (*climb) << endl;
		}
		if (finish.size() == 0) {
			throw new app_exception("[B] \'finish\' leafset is empty");
		}
		domainOfAttraction[finish]++;
//		DEBUG(cout << "Domain of attraction for split " << strmLeafset(finish) << " is increased to " << domainOfAttraction[finish] << endl); // XXXHERE
//		DEBUG(cout << "size of walks map = " << walks.size() << endl);
		return bestScore;
	}
	// Otherwise, we have an improving neighbour, so go to it and repeat:
//	if (b == current) {
//		return bestScore;
//	}
	DEBUG(cout << "Moving on to best neighbour: " << strmLeafset(b) << endl);
	bestScore = recursiveSteepestDescent(b, currentScore, finish, walkLength, os);
	Climb* climb = new Climb(currentSplit, getSplitError(currentSplit), finish, walkLength+1, bestScore);
	if (os.good()) {
		os << (*climb) << endl;
	}
	if (_verbose) {
		cout << (*climb) << endl;
	}
	++walkLength;
	walks[currentSplit] = climb;
	if (finish.size() == 0) {
		throw new app_exception("[C] \'finish\' leafset is empty");
	}
	domainOfAttraction[finish]++;
	DEBUG(cout << "size of walks map = " << walks.size() << endl);
	DEBUG(cout << "FINISHED recursiveSteepestDescent: " << strmLeafset(currentSplit) << endl);
	return bestScore;
}

void Project::setDataType(string dt) {
	/**
	 * So far this doesn't do anything -- but you never know; some time I might want to project
	 * characters on to R/Y coding or something.
	 */
	transform(dt.begin(), dt.end(), dt.begin(), ::tolower);
	if (!strcmp(dt.c_str(), "dna")) {
		cdType = cdt_DNA;
	} else if (!strcmp(dt.c_str(), "rna")) {
		cdType = cdt_RNA;
	} else if (!strcmp(dt.c_str(), "nucleotide")) {
		cdType = cdt_Nucleotide;
	} else if (!strcmp(dt.c_str(), "protein")) {
		cdType = cdt_Protein;
	} else if (!strcmp(dt.c_str(), "continuous")) {
		cdType = cdt_Continuous;
	}
}

void Project::setGraphColourParameters(double h, double s, double v) {
	graphColHue = h;
	graphColSaturation = s;
	graphColValue = v;
}

void Project::setGraphMaxScoreColour(double red, double green, double blue, double alpha) {
	graphColMaxScoreRed = red;
	graphColMaxScoreGreen = green;
	graphColMaxScoreBlue = blue;
	graphColMaxScoreAlpha = alpha;
//	cout << "Setting RGBA score for *maximum* score: " << red << ',' << green << ',' << blue << ',' << alpha << endl;
}

void Project::setGraphMinScoreColour(double red, double green, double blue, double alpha) {
	graphColMinScoreRed = red;
	graphColMinScoreGreen = green;
	graphColMinScoreBlue = blue;
	graphColMinScoreAlpha = alpha;
//	cout << "Setting RGBA score for *minimum* score: " << red << ',' << green << ',' << blue << ',' << alpha << endl;
}

void Project::run() {
	bool _debugging(true);
	if (splitsfile != "") {
		ifstream in(splitsfile);
			// this assumes a strict format: each split is a string of bits, then there's a comma, then a float.
		string str;
		char *tok;
		float val;
		string line;
		while (getline(in, line)) {
			char * cstr = new char [line.length()+1];
			std::strcpy (cstr, line.c_str());
			tok = strtok(cstr, ",");
			val = atof(strtok(NULL, ","));
			trueBranches[tok] = val;
			DEBUG(cout << "input split = " << tok << endl);
			delete[] cstr;
		};
		in.close();
//		for (auto spl : trueBranches) {
//			cout << spl.first << ":->" << spl.second << endl;
//		}
	}
	if (_doLandscapeAnalysis) {
		doLandscapeAnalysis();
	} else {
		calcPatternWeights();
		if (!_silent) {
			cout << "Calculating subflattening errors..." << endl;
		}
		calcAllSubflatteningErrors();
	}
	if (!_silent) {
		cout << "Done." << endl;
	}
}

void Project::setUpAlignmentPointers() {
	bool _debugging(true);
	DEBUG(cout << "A = " << A << endl);
	A = new string*[ntax];
	unsigned int i(0);
	for (auto taxonseq = rawAlignment.begin(); taxonseq != rawAlignment.end(); ++taxonseq) {
		A[i] = &(taxonseq->second);
		DEBUG(cout << "A[" << i << "] = " << A[i] << " = " << (*A[i]) << endl);
		++i;
	}
}


void Project::shortenSplit(leafset& a) {
	unsigned int card = count(a.begin(), a.end(), true);
	if (2*card > ntax) {
		a = flipLeafset(a);
	} else if ((2*card == ntax) && (a[0] == false)) {
		a = flipLeafset(a);	// for balanced splits, ensure the first (lexicographically) taxon is IN the split
	}
}

void Project::showSignedCounts() {
	for (pair<Pattern<unsigned int>, int> p : signedCountSums) {
		cout << p.first << "\t" << p.second << endl;
	}
}

void Project::showPatternCounts() {
	for (auto &p : patternCounts) {
		cout << p.first << "\t" << p.second << endl;
	}
	cout << "Total number of unique patterns found: " << patternCounts.size() << endl;
}

double Project::steepestDescentFromHere(leafset& current, unsigned int& walkLength) {
	if (walks.find(current) != walks.end()) {
		Climb* climb = walks.at(current);
		walkLength = climb->walkLength + 1;
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
//	bool _debugging(true);
	leafset current = getRandomSplit();
	unsigned int ignoreWalkLength;
	return steepestDescentFromHere(current, ignoreWalkLength);
}

string Project::strmSplitScoreMethod() {
	switch (sf) {
	case sf_cluster_distance:
		return "clust";
	case sf_singular_value:
		return "svd";
	case sf_parsimony:
		return "par";
	case sf_state_distribution:
		return "freq";
	default:
		return "";
	}
}


string Project::strmSplitAsTaxa(const leafset& split) {
	stringstream ss;
	ss << "{ ";
	auto t = taxon.begin();
	for (unsigned int i = 0; i < split.size(); ++i) {
		if (split[i]) {
			ss << *t << " ";
		}
		++t;
	}
	ss << "}";
	return ss.str();

}

void Project::putHalfSplitMask(vector<Pattern<char>>& idx, const flatbush::leafset& X) {

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
		for (unsigned int a = 1; a < numStates; ++a) {
			for (unsigned int left = 0; left < card; ++left) {
				// XXX check the value of both indices here: looks like an out of ranger error
				if (left >= position.size()) {
					cout << "left = " << left << "; position vector has size " << position.size() << endl;
					throw new app_exception("something's gone wrong with position vector.");
				}
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
	bool _debugging(true);
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

	// Get all the important input information from the file name:

	//

	captureAlignmentSettings(filename);
	DEBUG(cout << "setting project" << endl);
	parser->setProject(this);
	DEBUG(cout << "project set; now beginning parsing" << endl);
	parser->parse();
	DEBUG(cout << "checking sequence length... "; cout.flush());
	checkSequenceLength(true);
	DEBUG(cout << "sequence length = " << this->getSequenceLength() << "; ok\n");
	DEBUG(cout << "number of taxa = " << this->getNumTaxa() << "; ok\n");
//	DEBUG(cout << "compressing alignment... "; cout.flush());
//	compressAlignment();
	DEBUG(cout << "ok\n");
	checkSequenceLength(true);
	DEBUG(cout << "Parsing complete." << endl);
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
