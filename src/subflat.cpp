/*
// * subflat.cpp
 *
 *  Created on: 23 May 2016
 *      Author: mc7
 */
/**
 * TODO
 * TODO Indicate number of (potential) splits to check
 * TODO Output dimensions of each subflattening matrix?
 * TODO Give time started; estimate (if in verbose) how much longer to go...
 * TODO Perturb split and do loal search
 * TODO Check landscape of split space.
 * TODO Find out why there are these occasional issues with invalid input characters...
 * TODO Output (in level 3 diagnostic mode) the time taken for each split:
 * 	it looks like it's faster with more unbalanced splits, which is logical, but it would be good to measure it.
 * TODO Duh – I'm searching from each individidual split, so no duplicates there, but I'm not recording any other previous searches.
 * 	for each split encountered, register its final split found and how long it took to get there.
 */

#include <set>
#include <array>
#include <algorithm>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <map>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

//#include "../utility/mac.h"
#include "../utility/appexception.h"
#include "../utility/debugging.h"
#include "../utility/node.h"
#include "../utility/pattern.h"
#include "../utility/phylogeny.h"
#include "../utility/NEXUSParser.h"
#include "project.h"
#include "subflat.h"

/**
 * Alignment of sequences will be represented by pattern frequencies in a hash table.
 */

//const unsigned int k = 4;	// k is the size of the alphabet.

bool _debugging = false;
bool _verbose = false;
bool _silent = false;	// for silent running using batches of files or scripts.
using namespace std;
using namespace flatbush;

/*
 * Prototypes:
 */
unsigned int Fitch(const unsigned int* T, const leafset & C, unsigned int n);
std::string splash();
//void putHalfSplitMask(std::vector<Pattern<char>>& idx, leafset X);
std::vector<double> svdSplit(std::map<Pattern<unsigned int>, int> signedCount, unsigned int n, std::vector<bool> A);
//void putAllNontrivialHalfSplits(std::vector<leafset>& S, int ntax, bool _verbose);
void readNewickTree(flatbush::Node * n, const char* treestr);
void timeSVD();

string usage;

/**
 * I should use the SVDLIBC library for getting SVD values from sparse matrices as splitsup (Allman et al.) does,
 * if subflattenings are really any good.
 * See https://tedlab.mit.edu/~dr/SVDLIBC/ for download and installation.
 */
ofstream out;
Eigen::MatrixXi F;
Eigen::SparseMatrix<double> SF;

//void calcAllSubflatteningErrors(Project *proj, bool _verbose) {
//	unsigned int ntax = proj->getNumTaxa();
//	leafset right;
//	vector<Pattern<char>> rows;
//	vector<Pattern<char>> columns;
//	Pattern<char> pat(ntax, '.');
//	Pattern<unsigned int> idx;
//	Pattern<char> nullPattern(ntax, ' ');
//	vector<leafset>* splits = proj->getSplits();
//	for (leafset left : *splits) {
//		right = left;
//		right.flip();
//		putHalfSplitMask(rows, left);
//		putHalfSplitMask(columns, right);
//		if (_verbose) {
//			cout << "Creating sub-flattenings..." << endl;
//			cout.flush();
//		}
//		Eigen::MatrixXf F(rows.size(), columns.size());
//		if (_verbose) {
//			cout << "\t" << nullPattern << "\t";
//			for (auto col : columns) {
//				cout << setw(10 - ntax) << col << "\t";
//			}
//			cout << endl;
//		}
//
//	}
//
//}

void strmUsage() {
	stringstream ss;
	ss << "usage: flatbush [OPTIONS]\n";
	ss << "\t-f <filename>: input file name in fasta or nexus format.\n";
	ss << "\t\tNo default value --- without this, nothing much happens.\n";
	ss << "\t-l: do a landscape analysis.\n";
	ss << "\t-o <filename>: output file name for results.\n";
	ss << "\t\tDefault value is \"flatbush.out\".\n";
	ss << "\t-of <int>: output format for errors and splits.\n";
	ss << "\t\t0:\tbinary representation of each split, comma-separated-values (CSV);\n";
	ss << "\t\t1:\tbinary representation of each split, tab-delimited;\n";
	ss << "\t\t2:\tsplits as sets, comma-separated (CSV).\n";
	ss << "\t-s flag to turn on Silent mode for scripting:\n";
	ss << "\t\tDefault value is FALSE.\n";
	ss << "\t-v flag to turn on Verbose mode:\n";
	ss << "\t\tDefault value is FALSE.\n";
	usage = ss.str();
	cout << usage;
}

unsigned int Fitch(const unsigned int* T, const leafset & C, unsigned int n) {
/**
 * Calculate and return the parsimony length of leafset C, treated as a binary character, on tree T.
 */
	// lapsing into a little C here
	unsigned int* states = new (nothrow) unsigned int[2*n];	// this will be where the guts of the Fitch algorithm will go
	for (unsigned int i = 0; i < n; ++i) {
		states[i] = C[i] ? 2 : 1;
		cout << "states[" << i << "] = " << states[i] << endl;
	}
	unsigned int length = 0;
	for (unsigned int i = 0; i < 2*n-1; ++i) {
		if (T[i] > 2*n-1) {
			throw new app_exception("Tree index out of range");
		}
		cout << "states[T[" << i << "]] == " << states[T[i]] << endl;
		unsigned int a = (states[T[i]] == 0) ? states[i] : states[T[i]] & states[i];
//		cout << "\t<- " << states[T[i]] << endl;
		if (a == 0) {
			states[T[i]] |= states[i];
			cout << "state at parent node = " << states[T[i]] << endl;
			++length;
		} else {
			states[T[i]] = a;
		}
		cout << "states[T[" << i << "]] = " << states[T[i]] << endl;
	}
	/**
	 * Treating trees in the compact form of Joe F:
	 * [ 5 5 6 6 7 7 8 8 ]
	 * [ 0 1 2 3 4 5 6 7 ]
	 * corresponds to the tree
	 * 8---7---5---1
	 * |   |   +---0
	 * |   |
	 * |   +-------4
	 * |
	 * +---6-------3
	 *     +-------2
	 *
	 */
	delete [] states;
	return length;
}

/**
 * XXX Investigate fast calculation of just the first few SVs (mentioned in Allman, Kubatko & Rhodes)
 */
int main(int argn, char** argc) {
try {
	if (argn < 2) {
		strmUsage();
		return 1;
	}
//	DEBUG(_verbose = true);
//#ifndef DEBUGGING
//	_debugging = false;
//#endif
	/*
	 * if DEBUG flag is on, i.e., if it's been turned on by the preprocessor
	 * directive in the make command, then it makes sense to turn on _debugging
	 * *flag* by default.
	 */
	string fileName;
	string outFileName("flatbush.out");
	cout << "FlatBush City\n";
	Project proj;
	if (!_silent) {
		cout << "Reading input file..." << endl;
	}
	for (int i = 1; i < argn; ++i) {
		if (!strcmp(argc[i], "-v")) {
			_verbose = true;
		} else if (!strcmp(argc[i], "--debug")) {
			_debugging = true;
		} else if (!strcmp(argc[i], "-f")) {
			++i;
			fileName = argc[i];
			DEBUG(cout << "Read the input file" << endl);
			proj.read(fileName);
		} else if (!strcmp(argc[i], "-h")) {
			strmUsage();
		} else if (!strcmp(argc[i], "-l")) {
			proj.doLandscapeAnalysis();
			return 0;
		} else if (!strcmp(argc[i], "-o")) {
			outFileName = argc[i + 1];
		} else if (!strcmp(argc[i], "-of")) {
			++i;
			int n = atoi(argc[i]);
			if (n == fmt_binary_csv) {
				proj.outfmt = fmt_binary_csv;
			} else if (n== fmt_binary_tabbed) {
				proj.outfmt = fmt_binary_tabbed;
			} else {
				proj.outfmt = fmt_sets;
			}
		} else if (!strcmp(argc[i], "-r")) {
			proj.calcPatternWeights();
			proj.steepestDescentRandomStart();
			return 0;
		} else if (!strcmp(argc[i], "-s")) {
			++i;
			_silent = true;
		} else if (!strcmp(argc[i], "--seed")) {
			++i;
			srand48(atoi(argc[i]));
		} else if (!strcmp(argc[i], "--do-all-splits")) {
			++i;
			proj._doAllSplits = true;
		} else if (!strcmp(argc[i], "--buildtrees")) {
			++i;
			cout << "You've asked to build trees using this method, but it's not yet coded. Please be patient..." << endl;
		}
	}
	_silent = true;
	out.open(outFileName);
	if (out.good()) {
		cout << "Writing all output to file \'" << outFileName << "\'" << endl;
	} else if (!_verbose) {
		cout << "No output file and _verbose set to false. Please try again "
				<< "using -o <outfilename> to set the output file or -v to set verbose." << endl;
		return 1;
	}
	proj.calcPatternWeights();
	if (!_silent) {
		cout << "Calculating subflattening errors..." << endl;
	}
	proj.calcAllSubflatteningErrors();
	if (!_silent) {
		cout << "Done." << endl;
	}
	return 0;
} catch (exception* e) {
	cerr << "Exception: " << e->what() << endl;
}
}

string splash() {
	string str("FlatBush : a quirkily named program to calculate sub-flattenings.");
	str += "Build date: 20181210\nAuthor: M. A. Charleston\ne-mail: michael.charleston@utas.edu.au\n";
	return str;
}

//void putHalfSplitMask(vector<Pattern<char>>& idx, flatbush::leafset X) {
//
//	/**
//	 * put the "half" split patterns that can go into the row or column indices.
//	 * The format is like this:
//	 * 	..0..100..
//	 * which corresponds to the split (counting from 0, not 1) of { 2, 5, 6, 7 } | { 0, 1, 3, 4, 8, 9 }.
//	 * The pattern is the bit that changes:
//	 * 	..0..000..
//	 * 	..0..001..
//	 * 	..0..010..
//	 * 	etc.
//	 */
//	idx.clear();
//	Pattern<char> pat;
//	// all 0s:
//	unsigned int n = X.size();
//	for (unsigned int i = 0; i < n; ++i) {
//		char ch = '.';
//		if (X[i]) {
//			ch = '0';
//		}
//		pat.append(ch);
//	}
////	cout << pat << endl;
//	idx.push_back(pat);
////	unsigned int card = std::count_if(X.begin(), X.end(), bind1st(equal_to<bool>(), true));
//	unsigned int card = 0;
//	for (unsigned int i = 0; i < n; ++i) {
//		if (X[i]) {
//			++card;//std::count(X.begin(), X.end(), static_cast<bool>(true));
//		}
//	}
//	vector<unsigned int> position;
//	// in the example above, positions will be set to be the vector ( 2, 5, 6, 7 )
//	for (unsigned int i = 0; i < n; ++i) {
//		if (X[i]) {
//			position.push_back(i);
////			cout << "position += " << i << endl;
//		}
//	}
//	// singletons:
//	for (unsigned int a = 1; a < k; ++a) {
//		pat = '.';
//		for (unsigned int pos = 0; pos < card; ++pos) {
//			for (unsigned int left = 0; left < pos; ++left) {
//				pat[position[left]] = '0';
//			}
//			pat[position[pos]] = static_cast<char>('0' + a);
//			for (unsigned int right = pos+1; right < card; ++right) {
//				pat[position[right]] = '0';
//			}
////			cout << pat << endl;
//			idx.push_back(pat);
//		}
//	}
//}

/**
 *
 * This is the transformation from raw integer counts into the "Hadamard-transformed" counts:
 */
vector<double> svdSplit(map<Pattern<unsigned int>, int> signedCount, unsigned int n, vector<bool> A) {
	vector<double> diagonals;
//	Eigen::MatrixXf F;
	// construct the complement leaf set, B:
	vector<bool> B(n, false);
	for (unsigned int i = 0; i < n; ++i) {
		if (A.at(i)) {
			continue;	// i is in A
		}
		// i is not in A:
		B[i] = true;
	}
	// now create row and column indices:
	vector<Pattern<unsigned int>> colIdx;

	vector<Pattern<unsigned int>> rowIdx;

	return diagonals;	// don't mind a copy being made, not too onerous
}


void readNewickTree(Node * n, const char* treestr) {
	/**
	 * Read a tree in Newick format, as in "(A,(B,C));" and convert it to compact integer array.
	 * First number all the leaves:
	 * 	A -> 0
	 * 	B -> 1
	 * 	C -> 2
	 * then each pair of parentheses defines an internal node: number these ascending from the leaves:
	 * 	'(' ->
	 */

}

void timeSVD() {
	time_t start = clock();
	for (int i = 0; i < 1000; ++i) {
		Eigen::MatrixXf M = Eigen::MatrixXf::Random(20,30);
		Eigen::JacobiSVD<Eigen::MatrixXf> svd(M);
//		cout << "Rank of this randomly chosen matrix is " << lu.rank() << endl;
	}
	cout << "1000 matrix SVDs takes " << static_cast<float>(clock() - start) / CLOCKS_PER_SEC
			<< " seconds" << endl;
	Eigen::MatrixXf M = Eigen::MatrixXf::Random(20,30);
	Eigen::JacobiSVD<Eigen::MatrixXf> svd(M);
	cout << svd.singularValues() << endl;
}

