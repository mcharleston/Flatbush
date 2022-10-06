/*
 * project.h
 *
 *  Created on: 25 Jul 2016
 *      Author: mac
 */

#ifndef SRC_PROJECT_H_
#define SRC_PROJECT_H_

#include <set>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "../utility/pattern.h"
//#include "../utility/phylogeny.h"
#include "../utility/parser.h"
#include "subflat.h"
#include "Climb.h"
//#include "../utility/NEXUSParser.h"
//#include "../utility/CompactSequence.h"

namespace flatbush {

enum fileType {
	fastaFile,
	nexusFile,
	nullFile
};
enum CharacterDataType {
	cdt_Standard,
	cdt_DNA,
	cdt_RNA,
	cdt_Nucleotide,
	cdt_Protein,
	cdt_Continuous
};
enum DataType {
	alignmentData,
	patternWeightData
};
enum DistanceCorrection {
	dc_p,	// p distance
	dc_JC,	// Jukes-Cantor 69
	dc_HKY,	// Hasegawa-Kishino-Yano HKY85 correction
	dc_F81,
	dc_F84,
	dc_LogDet	// Log-Det method
};
enum OutputFormat {
	fmt_binary_csv,
	fmt_binary_tabbed,
	fmt_sets
};
enum splitFunction {
	sf_cluster_distance,
	sf_singular_value,
	sf_parsimony,
	sf_state_distribution
};

class Climb;
class Phylogeny;

//typedef std::map<Pattern<char>, unsigned int> PatternCounts;
typedef std::unordered_map<std::string, unsigned int> PatternCounts;
typedef std::unordered_map<std::string, double> PatternWeights;
typedef std::pair<leafset, leafset> edge;

size_t card(const leafset& l);

std::string strmLeafset(const leafset& vec);
std::string strmLeafsetBinary(const leafset& vec);

class Project {
private:
	std::string name;	/// name of this project; used as base name for file outputs, e.g. name-landscape.csv, name-graph.dot.
	std::string inputNotes;
	parsing::Parser* parser;
	fileType ft;
	DataType dType;
	CharacterDataType cdType;
	DistanceCorrection distanceMethod;
	std::vector<flatbush::leafset> splitsAsLeafsets;
	std::vector< std::set<std::string> > taxonSets;
	std::set<std::string> splitsAsStrings;

	// Data:
//	std::vector<Phylogeny> trees;
	size_t ntax;		/// number of taxa (= tips = leaves)
	std::map<std::string, int> leafLabels;
	std::unordered_map<std::string, float> trueBranches;	// both the splits and their branch lengths
	std::set<std::string> trueBSplits;
	float** D;	// distance matrix

	// Calculation parameters:
	splitFunction sf;

	// For alignments:
	std::map<std::string, std::string> rawAlignment;
	std::string** A;	// an array of pointers to strings; A[0] will contain the address of the sequence of the first (lex order) taxon
//	std::map<std::string, CompactSequence> compactAlignment;
	std::set<std::string> taxon;
	std::string substitutionModel; // if known; often provided as input from simulation software if it's in the filename as regex -([A-Z]+[0-9]+)
	float scaleFactor; // a useful input parameter, often provided as input from simulation software if it's in the filename as regex -[bh]([0-9]+\.[0-9]+)

	size_t seqLength;	/// in case of an alignment, the sequence length
	char gapChar = ' ';
	char missingChar = '?';
	char matchChar = '.';
	unsigned int numStates; // this is hard coded as 4 at present.
	std::string symbols;	// acceptable characters in alignment data, e.g. ACGT
	PatternCounts patternCounts;

	// For input pattern weights, e.g. probabilities
	PatternWeights patternWeights;

	// Calculation outputs:
	bool _hasPatternWeights;

	std::map<std::string, std::vector<int> > seqFreqs;
	std::map<Pattern<unsigned int>, int> signedCountSums;
	std::map<Pattern<unsigned int>, double> signedWeightSums;
	std::unordered_map<leafset, double> splitError;
//	std::unordered_map<leafset, unsigned int> parsimonySplitScore;
	std::unordered_map<leafset, Climb* > walks; // XXX HERE

	// Landscape features:
	unsigned int perts;
	unsigned int sampleSize;
	std::unordered_map<leafset, int> domainOfAttraction;	// XXX Need to have a good container here but don't want to pass objects, just pointers..
	std::vector<edge> climbEdges;
	std::string treefile; // name of an input tree file if known
	std::string splitsfile;	// name of an input file of splits if given

	// output fancy bits
	double graphColHue;
	double graphColSaturation;
	double graphColValue;

//	std::ostream os; // main output stream

public:
	bool _doAllSplits;
	bool _doLandscapeAnalysis;
	bool _suppress_labels;

	double graphColMaxScoreRed;
	double graphColMaxScoreGreen;
	double graphColMaxScoreBlue;
	double graphColMaxScoreAlpha;
	double graphColMinScoreRed;
	double graphColMinScoreGreen;
	double graphColMinScoreBlue;
	double graphColMinScoreAlpha;

	OutputFormat outfmt;
//	std::string landscapeFilename;
	static const std::map<char, int> encoding;

	virtual ~Project();
	Project();

	Project(std::string filename);

	void addKnownSplit(std::string s);
//	void addRawSequence(std::string&, std::string&);
	void addSequence(std::string, std::string);
	void addSplit(std::set<std::string> s);
	void addSplit(std::string s);
	void addTaxon(std::string id) { taxon.insert(id); rawAlignment[id]; }
	bool adjacent(const leafset& a, const leafset& b);
	PatternCounts bootstrapPatternCounts() const;
	void calcAllSubflatteningErrors();
	void calcDistanceMatrix(DistanceCorrection corr);
	unsigned int calcHammingDistance(const std::string* s, const std::string* t) const;
	float calcLogDetDistance(const std::string* s, const std::string* t) const;
	float calcMeanBetweenClusterDistance(const leafset& left) const;
	float calcMeanWithinClusterDistance(const leafset& left) const;
	void calcPatternWeights();
	double calcPDistance(const std::string* s, const std::string* t);
	void calcSeqFreqs();
	double calculateFlatteningMatrixError(leafset split);
	void captureAlignmentSettings(std::string str);
	void checkSequenceLength(bool nothrow=false);
	leafset chooseRandomSplit();
	void combine(Pattern<char>& result, const Pattern<char> x, const Pattern<char> y);
	void compressAlignment();
	void countAllSignedSums(bool _verbose);
	unsigned int countMatchingSplits(std::unordered_map<std::string, float>, std::set<std::string>);
	void countPatterns();
	unsigned int getCard(const leafset& leaves);
	int countSignedSum(Pattern<unsigned int> pat);
	std::string describeSplit(leafset pat);
	void doLandscapeAnalysis();
	void exportLandscape();
	double findBestNeighbouringSplit(leafset& bestNeighbour, leafset L, double score);
	leafset flipLeafset(const leafset& L);
	void gatherAllNontrivialHalfSplits(bool _verbose);
	void gatherTaxonLabels();
//	Alignment* getAlignment() { return &A; }
	char getAlignmentMatchChar() { return matchChar; }
	DataType getDataType() const { return dType; }
	unsigned int getNumStates() const { return numStates; }
	size_t getNumTaxa() const { return ntax; }
	PatternCounts& getPatternCounts() { return patternCounts; }
	PatternWeights* getPatternWeights() { return &patternWeights; }
	leafset getRandomSplit();
	std::map<std::string, std::string>* getRawAlignment() { return &rawAlignment; }
	std::string& getRawSequence(const std::string& str) { return rawAlignment.at(str); }
	std::string& getSequence(unsigned int seqID);
	const std::string& getSequence(unsigned int seqID) const;
	std::string& getSequence(const std::string& seqID);
	const std::string& getSequence(const std::string& seqID) const;
	unsigned int getSequenceLength() const { return seqLength; }
	std::map<std::string, std::string>* getSequenceMap() { return &rawAlignment; }
	std::map<Pattern<unsigned int>, int>& getSignedSums() { return signedCountSums; }
	double getSplitError(const flatbush::leafset& split);
	double getSplitErrorClustering(const flatbush::leafset& left);
	double getSplitErrorStateDistribution(const flatbush::leafset& split);
	double getSplitErrorSVD(const flatbush::leafset& split);
	double getSplitErrorParsimony(const flatbush::leafset& split);
//	double getSplitErrorStateDistribution(const flatbush::leafset& split);
	std::vector<leafset>* getSplits();
	std::string& getStates() { return symbols; }
	unsigned int HammingDistance(const leafset& a, const leafset &b);
	bool hasPatternWeights() const { return _hasPatternWeights; }
	bool hasSequence(std::string seqID);
	void initialise();
	bool isTrivialSplit(const leafset& split);
	void putAllAdjacentSplits(std::vector<flatbush::leafset>& N, const flatbush::leafset L);
	void putAllHalfSplitPatterns(std::vector<Pattern<char>>& idx, leafset X);
	void putHalfSplitMask(std::vector<Pattern<char>>& idx, const leafset& X);
	void read(const std::string& fileName);
	double recursiveSteepestDescent(const leafset& current, double currentScore, leafset& finish,
			unsigned int& walkLength, std::ostream& os);
	void run();
//	void setAlignmentDataType(std::string& dt) { A.setDataType(dt); }
	void setAlignmentGapChar(char ch) { gapChar = ch; }
	void setAlignmentMatchChar(char ch) { matchChar = ch; }
	void setAlignmentMissingChar(char ch) { missingChar = ch; }
//	void setDataType(const std::string& dt);
//	void setDataType(const std::string& s);
	void setDataType(std::string dt);
	void setGapChar(char c) { gapChar = c; }
	void setGraphColourParameters(double red, double green, double blue);
	void setGraphMaxScoreColour(double red, double green, double blue, double alpha);
	void setGraphMinScoreColour(double red, double green, double blue, double alpha);
	void setMissingChar(char c) { missingChar = c; };
	void setName(const std::string& base) { name = base; }
	void setNTax(int n) { ntax = n; }
	void setNumSites(int c) { seqLength = c; }
	void setNumTaxa(int n) { ntax = n; }
	void setPerturbations(unsigned int p) { perts = p; }
	void setSampleSize(int n) { sampleSize = n; }
	void setScoreMethod(splitFunction split_function) { sf = split_function; }
	void setSeqLength(int c) { seqLength = c; }
	void setSplitsFile(const std::string& splfile) { splitsfile = splfile; }
	void setSymbols(const std::string& sym) { symbols = sym; }
	void setTreeFile(const std::string& trfile) { treefile = trfile; }
	void setUpAlignmentPointers();
	void showPatternCounts();
	void showSignedCounts();
	void shortenSplit(leafset &a);
	double steepestDescentFromHere(leafset& start, unsigned int& walkLength);
	double steepestDescentRandomStart();
	std::string strmSplitScoreMethod();
	std::string strmSplitAsTaxa(const leafset& split);
};

void combine(Pattern<char>& result, const Pattern<char> x, const Pattern<char> y);
Pattern<unsigned int>& encode(Pattern<unsigned int>& result, const Pattern<char>& pat);

} /* namespace flatbush */

#endif /* SRC_PROJECT_H_ */
