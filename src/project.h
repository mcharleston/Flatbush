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
#include "../utility/pattern.h"
//#include "../utility/phylogeny.h"
#include "../utility/parser.h"
#include "subflat.h"
#include "Climb.h"
//#include "../utility/NEXUSParser.h"
#include "../utility/CompactSequence.h"

namespace flatbush {

enum fileType {
	fastaFile,
	nexusFile,
	nullFile
};
enum DataType {
	alignmentData,
	patternWeightData
};
enum OutputFormat {
	fmt_binary_csv,
	fmt_binary_tabbed,
	fmt_sets
};

class Phylogeny;
class Climb;

//typedef std::map<Pattern<char>, unsigned int> PatternCounts;
typedef std::unordered_map<CompactSequence, unsigned int> PatternCounts;
typedef std::unordered_map<CompactSequence, double> PatternWeights;

size_t card(const leafset& l);
std::string strmLeafset(const leafset& vec);
std::string strmLeafsetBinary(const leafset& vec);

class Project {
private:
	std::string name;	/// name of this project
	parsing::Parser* parser;
	fileType ft;
	DataType dType;
	std::vector<flatbush::leafset> splits;
	std::vector< std::set<std::string> > taxonSets;
	std::set<std::string> knownSplits;

	// Data:
//	std::vector<Phylogeny> trees;
	size_t ntax;		/// number of taxa (= tips = leaves)
	std::map<std::string, int> leafLabels;

	// For alignments:
	std::map<std::string, std::string> rawAlignment;
	std::map<std::string, CompactSequence> compactAlignment;
	size_t seqLength;	/// in case of an alignment, the sequence length
	char gapChar;
	char missingChar;
	unsigned int numStates; // this is hard coded as 4 at present.
	PatternCounts patternCounts;

	// For input pattern weights, e.g. probabilities
	PatternWeights patternWeights;

	// Calculation outputs:
	bool _hasPatternWeights;

	std::map<Pattern<unsigned int>, int> signedCountSums;
	std::map<Pattern<unsigned int>, double> signedWeightSums;
	std::unordered_map<leafset, double> splitError;
	std::unordered_map<leafset, Climb> walks;
public:
	OutputFormat outfmt;
	bool _doAllSplits;
	static const std::map<char, int> encoding;

	virtual ~Project();
	Project();

	Project(std::string filename);

	void addKnownSplit(std::string s);
	void addRawSequence(std::string, std::string);
	void addSequence(std::string, std::string);
	void addSplit(std::set<std::string> s);
	void addSplit(std::string s);
	void addTaxon(std::string id) { compactAlignment[id]; }
	bool adjacent(const leafset& a, const leafset& b);
	PatternCounts bootstrapPatternCounts() const;
	void calcAllSubflatteningErrors();
	void calcPatternWeights();
	double calculateFlatteningMatrixError(leafset split);
	void checkSequenceLength(bool nothrow=false);
	void compressAlignment();
	void countAllSignedSums(bool _verbose);
	void countPatterns();
	int countSignedSum(Pattern<unsigned int> pat);
	std::string describeSplit(leafset pat);
	void doLandscapeAnalysis();
	void exportLandscape(std::string & filename);
	double findBestNeighbouringSplit(leafset& bestNeighbour, leafset L, double score);
	void gatherAllNontrivialHalfSplits(bool _verbose);
	void gatherTaxonLabels();
//	Alignment* getAlignment() { return &A; }
	DataType getDataType() const { return dType; }
	unsigned int getNumStates() const { return numStates; }
	size_t getNumTaxa() const { return ntax; }
	PatternCounts& getPatternCounts() { return patternCounts; }
	PatternWeights* getPatternWeights() { return &patternWeights; }
	leafset getRandomSplit();
	std::map<std::string, std::string>* getRawAlignment() { return &rawAlignment; }
	std::string& getRawSequence(const std::string& str) { return rawAlignment.at(str); }
	CompactSequence getSequence(const std::string& seqID);
	unsigned int getSequenceLength() const { return seqLength; }
	std::map<std::string, CompactSequence>* getSequenceMap() { return &compactAlignment; }
	std::map<Pattern<unsigned int>, int>& getSignedSums() { return signedCountSums; }
	double getSplitError(const flatbush::leafset& split);
	std::vector<leafset>* getSplits();
	unsigned int HammingDistance(const leafset& a, const leafset &b);
	bool hasPatternWeights() const { return _hasPatternWeights; }
	bool hasSequence(std::string seqID);
	bool isTrivialSplit(const leafset& split);
	void putAllAdjacentSplits(std::vector<flatbush::leafset>& N, const flatbush::leafset L);
	void putAllHalfSplitPatterns(std::vector<Pattern<char>>& idx, leafset X);
	void putHalfSplitMask(std::vector<Pattern<char>>& idx, leafset X);
	void read(const std::string& fileName);
	double recursiveSteepestDescent(const leafset& current, double currentScore, leafset& finish, unsigned int& walkLength);
//	void setAlignmentDataType(std::string& dt) { A.setDataType(dt); }
	void setAlignmentGapChar(char ch) { gapChar = ch; }
	void setAlignmentMissingChar(char ch) { missingChar = ch; }
//	void setDataType(const std::string& dt);
//	void setDataType(const std::string& s);
	void setDataType(DataType t) { dType = t; }
	void setGapChar(char c) { gapChar = c; }
	void setMissingChar(char c) { missingChar = c; };
	void setNTax(int n) { ntax = n; }
	void setNumSites(int c) { seqLength = c; }
	void setNumTaxa(int n) { ntax = n; }
	void setSeqLength(int c) { seqLength = c; }
	void showPatternCounts();
	void showSignedCounts();
	double steepestDescentFromHere(leafset& start, unsigned int& walkLength);
	double steepestDescentRandomStart();
};

void combine(Pattern<char>& result, const Pattern<char> x, const Pattern<char> y);
Pattern<unsigned int>& encode(Pattern<unsigned int>& result, const Pattern<char>& pat);

} /* namespace flatbush */

#endif /* SRC_PROJECT_H_ */
