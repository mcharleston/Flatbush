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
#include <chrono>
#include <ctime>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <map>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "../utility/mac.h"
#include "../utility/tools.h"
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
const unsigned int kNumStates = 4;

bool _debugging = false;
bool _save_climbs(false);
bool _silent(false);	// for silent running using batches of files or scripts.
bool _show_all_details(false);
bool _show_taxon_names(true);	// if possible print taxon names in full in splits
bool _verbose(false);
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

unsigned short NucsAsBits[256];

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


// Nucleotide ambiguity codes from http://insilico.ehu.es/restriction/Nucleotide_ambiguity_code.html
void initialise() {
	// non-ambiguous
	NucsAsBits['a'] = 1;
	NucsAsBits['c'] = 2;
	NucsAsBits['g'] = 4;
	NucsAsBits['t'] = 8;
	NucsAsBits['A'] = 1;
	NucsAsBits['C'] = 2;
	NucsAsBits['G'] = 4;
	NucsAsBits['T'] = 8;
	// ambiguous:
	NucsAsBits['y'] = 10;	// Y = Pyrimidine: complement R
	NucsAsBits['Y'] = 10;
	NucsAsBits['r'] = 5;		// R = Purine: complement T
	NucsAsBits['R'] = 5;
	NucsAsBits['w'] = 9;		// W = Weak: complement S
	NucsAsBits['W'] = 9;
	NucsAsBits['s'] = 6;		// S = Strong: complement W
	NucsAsBits['S'] = 6;
	NucsAsBits['k'] = 12;	// K = keto: complement M
	NucsAsBits['K'] = 12;
	NucsAsBits['m'] = 3;		// M = amino: complement K
	NucsAsBits['M'] = 3;
	NucsAsBits['d'] = 13;	// D = not C: complement C
	NucsAsBits['D'] = 13;
	NucsAsBits['v'] = 7;		// V = not T: complement T
	NucsAsBits['V'] = 7;
	NucsAsBits['h'] = 11;	// H = not G: complement G
	NucsAsBits['H'] = 11;
	NucsAsBits['b'] = 14;	// B = not A: complement A
	NucsAsBits['B'] = 14;
	NucsAsBits['n'] = 15;	// N/X unknown
	NucsAsBits['N'] = 15;
	NucsAsBits['x'] = 15;
	NucsAsBits['X'] = 15;
}

unsigned int numStatesFromNucCode[16] = {
		[0] = 0,		// 0000
		[1] = 1,		// 0001
		[2] = 1,		// 0010
		[3] = 2,		// 0011
		[4] = 1,		// 0100
		[5] = 2,		// 0101
		[6] = 2,		// 0110
		[7] = 3,		// 0111
		[8] = 1,		// 1000
		[9] = 2,		// 1001
		[10] = 2,	// 1010
		[11] = 3,	// 1011
		[12] = 2,	// 1100
		[13] = 3,	// 1101
		[14] = 3,	// 1110
		[15] = 4		// 1111
};

void strmUsage() {
	stringstream ss;
	ss << "usage: flatbush [OPTIONS]\n";
	ss << "\t-h:\t stream this usage and quit\n";
	// -i input file
	ss << "\t-i <filename>: input file name in Fasta or Nexus format.\n";
	ss << "\t\tNo default value --- without this, nothing much happens.\n";
	// -l landscape
	ss << "\t-l: do a landscape analysis: from each split or from randomly selected splits, perform a steepest descent.\n";
	ss << "\t\tCreate an output graph in dot format with all the descents.\n";
	// -m scoring method
	ss << "\t-m <method>: scoring method for splits (default is subflattening SVD score), where <method> is one of\n";
	ss << "\t\tsvd: Singular Value Decomposition score of subflattening matrix;\n";
	ss << "\t\tparsimony|par: par - parsimony; summed minimum parsimony score of each site pattern on any tree containing the split\n";
	ss <<	"\t\tclust - cluster distance; sum of mean within-part Jukes-Cantor distance for each side of the split,\n";
	ss << "\t\t\tless the mean distance between taxon between sides of the split\n";
	ss << "\t\tfre - distance based on Euclidean distance between vectors of mean base frequencies on either side of the split\n";
	ss << "\t\tDefault method: par.\n";
	ss << "\t-n <string>: project base name from which all file names will be derived\n";
	// -o output project base file name
	ss << "\t-o <filename>: output file name for results.\n";
	ss << "\t\tDefault value is \"flatbush.out\".\n";
	// -of output format
	ss << "\t-of <int>: output format for errors and splits.\n";
	ss << "\t\t0:\tbinary representation of each split, comma-separated-values (CSV);\n";
	ss << "\t\t1:\tbinary representation of each split, tab-delimited;\n";
	ss << "\t\t2:\tsplits as sets, comma-separated (CSV).\n";
	ss << "\t\tDefault: 0 (binary, csv)\n";
	// -s silent mode
	ss << "\t-p[ert]: perturbation type(s) to use:\n";
	ss << "\t\t1 for single-taxon move only; 2 for swap only (retaining split sizes); 3 for both.\n";
	ss << "\t\tDefault: 2\n";
	ss << "\t-ss <int>: set the sample size to <int> value\n";
	ss << "\t\tDefault: 1000\n";
	// -v verbose
	ss << "\t-v flag to turn on Verbose mode:\n";
	ss << "\t\tDefault: off\n";
//	ss << "\t--hsv <float> <float> <float>:\tinput HSV (hue, saturation, value) colour parameters\n";
	ss << "\t-rgbamax <float> <float> <float> <float>: red, green, blue, alpha values respectively for the *minimum* (WORST) score\n";
	ss << "\t\tDefaults: 0.0 0.0 0.5 0.25\n";
	ss << "\t-rgbamin <float> <float> <float> <float>: red, green, blue, alpha values respectively for the *maximum* (BEST) score\n";
	ss << "\t\tDefaults: 0.0 1.0 0.0 0.0\n";
	ss << "\t-rgbmax <float> <float> <float>: red, green, blue values respectively for the *maximum* (WORST) score\n";
	ss << "\t\tDefaults: 0.0 1.0 0.0 (with alpha = 0.0)\n";
	ss << "\t-rgbmin <float> <float> <float>: red, green, blue values respectively for the *minimum* (BEST) score\n";
	ss << "\t\tDefaults: 0.0 0.0 0.5 (with alpha = 0.25)\n";
	ss << "\t--buildtrees: build trees via divide-and-conquer (not yet implemented!)\n";
	ss << "\t--debug: if set and FlatBush compiled with -DDEBUGGING defined, turn on debugging messages for problem diagnosis\n";
	ss << "\t\tDefault: off\n";
	ss << "\t--do-all-splits: systematically steepest descent from each split. If ntax>16 this will be over-ridden and splits will be sampled.\n";
	ss << "\t\tDefault: off\n";
	ss << "\t--save-climbs: create a csv file <project name>-climbs.csv with starting split & its score, final split and its score,\n";
	ss << "\t\tsize of the split (cardinality of smaller part) and path length. Potentially a big file.\n";
	ss << "\t--sdrs: do a single Steepest Descent Random Start (useful for debugging but not much else)\n";
	ss << "\t--seed <int>: set the random number seed, for repeatability\n";
	ss << "\t\tDefault: seed set from the system clock\n";
	ss << "\t--show-all-details: show all the split errors:\n";
	ss << "\t\tthis will potentially take a much longer time and fill up your terminal with junk you didn't really want.\n";
	ss << "\t--silent flag to turn on Silent mode, useful for bash scripting\n";
	ss << "\t\tDefault: off\n";
	ss << "\t--splitsfile <filename>: load a set of splits in as \"true\" splits\n";
	ss << "\t\tformat: a list of splits as binary strings with a floating-point number for each, such as branch length.\n";
	ss << "\t\te.g. \"00111, 0.5\" to define split {0,1} | {2,3,4} with value 0.5\n";
	ss << "\t--[no-]taxon-names: --taxon-names if possible, print splits with taxon names;\n";
	ss << "\t                    --no-taxon-names: do not use taxon names in splits.\n";
	ss << "\t\tDefault: off\n";
	ss << "\t--suppress-node-labels: instead of showing splits in the landscape graph, just show dots.\n";
	ss << "\t\tDefault: off\n";
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

typedef std::chrono::milliseconds millisecond;
typedef std::chrono::duration<float> fsec;

/**
 * XXX Investigate fast calculation of just the first few SVs (mentioned in Allman, Kubatko & Rhodes)
 */
int main(int argn, char** argc) {
try {
//	DEBUG(_verbose = true);
//#ifndef DEBUGGING
//	_debugging = false;
//#endif
	/*
	 * if DEBUG flag is on, i.e., if it's been turned on by the preprocessor
	 * directive in the make command, then it makes sense to turn on _debugging
	 * *flag* by default.
	 */
	auto wallclock_start = std::chrono::system_clock::now();	// "actual" time for noting when the program is run
	auto steady_start = std::chrono::steady_clock::now();	// "steady" clock won't ever be adjusted down so ideal for accurate calculating duration of run
	std::time_t start_time = std::chrono::system_clock::to_time_t(wallclock_start);

	if (argn < 2) {
		strmUsage();
		return (0);
	}
	bool _debugging(true);
	string fileName;
	string outFileName("flatbush.out");

	bool _testingHex(false);
	if (_testingHex) {
		double r, a;
		for (r = 0.0; r <= 1.0; r = r+0.1) {
			for (a = 0.0; a <= 1.0; a += 0.1) {
				cout << "r,a = " << r << ',' << a << " -> " << unitTo256Hex(r) << unitTo256Hex(a) << endl;
			}
		}
		return(0);
	}
	initialise();

	cout << "FlatBush";
	for (int i = 1; i < argn; ++i) {
		cout << ' ' << argc[i];
	}
	cout << ", start=" << strtok(std::ctime(&start_time), "\n") << ", ";
	Project proj;
	proj.setGraphColourParameters(0.4,0.8,0.3);
//	if (!_silent) {
//		cout << "Reading input file..." << endl;
//	}
	srand48((long int)time(NULL));

	// **GET ALL THE ARGUMENTS FIRST**
	for (int i = 1; i < argn; ++i) {
		if (!strcmp(argc[i], "-v")) {
			_verbose = true;
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
		} else if (!strcmp(argc[i], "-h")) {
			strmUsage();
			return(0);
		} else if (!strcmp(argc[i], "-i")) {
			++i;
			fileName = argc[i];
			DEBUG(cout << "Read the input file" << endl);
			proj.read(fileName);
		} else if (!strcmp(argc[i], "-l")) {
			proj._doLandscapeAnalysis = true;
		} else if (!strcmp(argc[i], "-m")) {
			++i;
			if (!strncmp(argc[i], "par", 3)) {
				// just test first 3 letters of "parsimony"
				DEBUG(cout << "Using Parsimony score for each split\n");
				proj.setScoreMethod(sf_parsimony);
			} else if (!strncmp(argc[i], "clust", 3)) {
				DEBUG(cout << "Using Cluster Distance score for each split\n");
				proj.setScoreMethod(sf_cluster_distance);
			} else if (!strncmp(argc[i], "svd", 3)) {
				DEBUG(cout << "Using Singular Value Decomposition score for each split\n");
				proj.setScoreMethod(sf_singular_value);
			} else if (!strncmp(argc[i], "fre", 3)) {
				DEBUG(cout << "Using base frequency distribution comparison for each split\n");
				proj.setScoreMethod(sf_state_distribution);
			}
		} else if (!strcmp(argc[i], "-n")) {
			++i;
			proj.setName(argc[i]); // over-rides the project name from the default
		} else if (!strcmp(argc[i], "-pert") || !strcmp(argc[i], "-p")) {
			++i;
			proj.setPerturbations(atoi(argc[i]));	// 1 for single-taxon move only; 2 for swap only; 3 for both.
//		} else if (!strcmp(argc[i], "-hsv")) { // hue/saturation/value colours for graph output
//			++i;
//			double h = atof(argc[i]);
//			++i;
//			double s = atof(argc[i]);
//			++i;
//			double v = atof(argc[i]);
//			proj.setGraphColourParameters(h, s, v);
//			DEBUG(cout << "colour parameters are HSV=(" << h << ',' << s << ',' << v << ")" << endl);
		} else if (!strcmp(argc[i], "--silent")) { // silent mode
			_silent = true;
		} else if (!strcmp(argc[i], "-ss")) { // sample size
			++i;
			proj.setSampleSize(atoi(argc[i]));
			DEBUG(cout << "Sample size = " << atoi(argc[i]) << endl);
		} else if (!strcmp(argc[i], "--buildtrees")) {
			cout << "You've asked to build trees using this method, but it's not yet coded. Please be patient..." << endl;
//		} else if (!strcmp(argc[i], "--comptrees")) {
//			++i;
//			while (!strcmp(argc[i+1], ";")) {
//				++i;
//				parsing::NEXUSParser NP();
//				Phylogeny* T(nullptr);
//
//			}
		} else if (!strcmp(argc[i], "--debug")) {
			_debugging = true;
		} else if (!strcmp(argc[i], "--do-all-splits")) {
			proj._doAllSplits = true;
		} else if (!strcmp(argc[i], "--hsv")) {
			++i;
			double h = atof(argc[i]);
			++i;
			double s = atof(argc[i]);
			++i;
			double v = atof(argc[i]);
			proj.setGraphColourParameters(h, s, v);
			DEBUG(cout << "colour parameters are HSV=(" << h << ',' << s << ',' << v << ")" << endl);
		} else if (!strcmp(argc[i], "-rgbamax")) {
			++i;
			double r = atof(argc[i]);
			++i;
			double g = atof(argc[i]);
			++i;
			double b = atof(argc[i]);
			++i;
			double a = atof(argc[i]);
			proj.setGraphMaxScoreColour(r, g, b, a);
		} else if (!strcmp(argc[i], "-rgbamin")) {
			++i;
			double r = atof(argc[i]);
			++i;
			double g = atof(argc[i]);
			++i;
			double b = atof(argc[i]);
			++i;
			double a = atof(argc[i]);
			proj.setGraphMinScoreColour(r, g, b, a);
		} else if (!strcmp(argc[i], "-rgbmax")) {
			++i;
			double r = atof(argc[i]);
			++i;
			double g = atof(argc[i]);
			++i;
			double b = atof(argc[i]);
			++i;
			proj.setGraphMaxScoreColour(r, g, b, proj.graphColMaxScoreAlpha);
		} else if (!strcmp(argc[i], "-rgbmin")) {
			++i;
			double r = atof(argc[i]);
			++i;
			double g = atof(argc[i]);
			++i;
			double b = atof(argc[i]);
			++i;
			proj.setGraphMinScoreColour(r, g, b, proj.graphColMinScoreAlpha);
		} else if (!strcmp(argc[i], "--no-taxon-names")) {
			_show_taxon_names = false;
		} else if (!strcmp(argc[i], "--save-climbs")) {
			_save_climbs = true;
		} else if (!strcmp(argc[i], "--seed")) {
			++i;
			srand48(atoi(argc[i]));
		} else if (!strcmp(argc[i], "--sdrs")) {
			proj.calcPatternWeights();
			proj.steepestDescentRandomStart();	// this stops execution after a single random start: mainly for debugging.
			return 0;
		} else if (!strcmp(argc[i], "--show-all-details")) {
			_show_all_details = true;
			cout << "You have opted to show all the split errors: ";
			cout << "this will potentially take a much longer time and fill up your terminal with crap you didn't really want.";
		} else if (!strcmp(argc[i], "--splitsfile")) {
			++i;
			proj.setSplitsFile(argc[i]);
		} else if (!strcmp(argc[i], "--suppress-node-labels")) {
			proj._suppress_labels = true;
		} else if (!strcmp(argc[i], "--taxon-names")) {
			_show_taxon_names = true;
		} else {
			cout << "I can't read this argument: please try again, bearing in mind what I expect:\n";
			cout << "Argument supplied: " << argc[i] << endl;
			strmUsage();
			throw new app_exception("Input argument parsing error.");
			return(1);
		}
	}
	_silent = true;
	proj.run();

	auto steady_end = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration<float>(steady_end - steady_start);
//	millisecond duration = std::chrono::duration_cast<millisecond>(elapsed);
	cout << endl << "Elapsed_time(s)=" << elapsed.count() << endl;
	return 0;
} catch (exception* e) {
	cerr << "Exception: " << e->what() << endl;
}
}

string splash() {
	string str("FlatBush : a quirkily named program to calculate sub-flattenings.");
	str += "Build date: 20200929\nAuthor: M. A. Charleston\ne-mail: michael.charleston@utas.edu.au\n";
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

