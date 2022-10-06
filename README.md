# FlatBush
A C++ program that calculates "subflattening" values to find promising bipartitions of phylogenetic data.

## Installation
Download the code from src/ and utility/ and compile with a C++11- compliant compiler such as gcc, e.g.:
```
> g++ -std=c++11 src/*.cpp utility/*.cpp -O3 -o FlatBush
```

Alternatively to build FlatBush in "DEBUG" mode, 
```
> g++ -std=c++11 src/*.cpp utility/*.cpp -DDEBUGGING -o FlatBush
```

## Running *FlatBush*
FlatBush is run from the command-line with optional arguments, but one obviously needed argument is a set of aligned sequence data.
Running FlatBush with no arguments
```
> ./FlatBush
```
yields
```
FlatBush -h, start=Thu Oct  6 16:51:24 2022, usage: flatbush [OPTIONS]
	-h:	 stream this usage and quit
	-i <filename>: input file name in Fasta or Nexus format.
		No default value --- without this, nothing much happens.
	-l: do a landscape analysis: from each split or from randomly selected splits, perform a steepest descent.
		Create an output graph in dot format with all the descents.
	-m <method>: scoring method for splits (default is subflattening SVD score), where <method> is one of
		svd: Singular Value Decomposition score of subflattening matrix;
		parsimony|par: par - parsimony; summed minimum parsimony score of each site pattern on any tree containing the split
		clust - cluster distance; sum of mean within-part Jukes-Cantor distance for each side of the split,
			less the mean distance between taxon between sides of the split
		fre - distance based on Euclidean distance between vectors of mean base frequencies on either side of the split
		Default method: par.
	-n <string>: project base name from which all file names will be derived
	-o <filename>: output file name for results.
		Default value is "flatbush.out".
	-of <int>: output format for errors and splits.
		0:	binary representation of each split, comma-separated-values (CSV);
		1:	binary representation of each split, tab-delimited;
		2:	splits as sets, comma-separated (CSV).
		Default: 0 (binary, csv)
	-p[ert]: perturbation type(s) to use:
		1 for single-taxon move only; 2 for swap only (retaining split sizes); 3 for both.
		Default: 2
	-ss <int>: set the sample size to <int> value
		Default: 1000
	-v flag to turn on Verbose mode:
		Default: off
	-rgbamax <float> <float> <float> <float>: red, green, blue, alpha values respectively for the *minimum* (WORST) score
		Defaults: 0.0 0.0 0.5 0.25
	-rgbamin <float> <float> <float> <float>: red, green, blue, alpha values respectively for the *maximum* (BEST) score
		Defaults: 0.0 1.0 0.0 0.0
	-rgbmax <float> <float> <float>: red, green, blue values respectively for the *maximum* (WORST) score
		Defaults: 0.0 1.0 0.0 (with alpha = 0.0)
	-rgbmin <float> <float> <float>: red, green, blue values respectively for the *minimum* (BEST) score
		Defaults: 0.0 0.0 0.5 (with alpha = 0.25)
	--buildtrees: build trees via divide-and-conquer (not yet implemented!)
	--debug: if set and FlatBush compiled with -DDEBUGGING defined, turn on debugging messages for problem diagnosis
		Default: off
	--do-all-splits: systematically steepest descent from each split. If ntax>16 this will be over-ridden and splits will be sampled.
		Default: off
	--save-climbs: create a csv file <project name>-climbs.csv with starting split & its score, final split and its score,
		size of the split (cardinality of smaller part) and path length. Potentially a big file.
	--sdrs: do a single Steepest Descent Random Start (useful for debugging but not much else)
	--seed <int>: set the random number seed, for repeatability
		Default: seed set from the system clock
	--show-all-details: show all the split errors:
		this will potentially take a much longer time and fill up your terminal with junk you didn't really want.
	--silent flag to turn on Silent mode, useful for bash scripting
		Default: off
	--splitsfile <filename>: load a set of splits in as "true" splits
		format: a list of splits as binary strings with a floating-point number for each, such as branch length.
		e.g. "00111, 0.5" to define split {0,1} | {2,3,4} with value 0.5
	--[no-]taxon-names: --taxon-names if possible, print splits with taxon names;
	                    --no-taxon-names: do not use taxon names in splits.
		Default: off
	--suppress-node-labels: instead of showing splits in the landscape graph, just show dots.
		Default: off
```

A slightly more involved example is
```
> ./FlatBush -i ../filo10Taxon.100.fst -n sim10 -l; xdot sim10-graph.dot &
```
... which 

* invokes FlatBush 
* to operate on the FASTA-format file filo10Taxon.100.fst in a directory above the current location of FlatBush, 
* setting the prioject name to "sim10",
* doing a landscape analysis (and for this few taxa, that means checking all splits),
* and then invoking the "xdot" Linux command to open the resulting graph for display.

### Another example

A slightly more involved example is the following, using input FASTA format file n10-1-1.fst and splits file n10-1.bsplits to calculate the landscape under the "cluster distance" score:
```
> ./FlatBush -i n10-1-1.fst --suppress-node-labels -n n10 --save-climbs -m clust -rgbamax 0.2 1.0 0.3 0.5 -rgbamin 1.0 0 1 0.8 --splitsfile n10-1.bsplits -l
```

The output on the command-line is
```
FlatBush -i n10-1-1.fst --suppress-node-labels -n n10 --save-climbs -m clust -rgbamax 0.2 1.0 0.3 0.5 -rgbamin 1.0 0 1 0.8 --splitsfile n10-1.bsplits -l, start=Thu Oct  6 17:09:04 2022, 
Elapsed_time(s)=0.0431931
```
and which outputs the files 

#### n10-climbs.csv
beginning
```
StartSplit,StartScore,StopSplit,StopScore,SplitSize,PathLength
{ 0 1 4 7 8 },-0.0285387,{ 0 1 4 7 8 },-0.0285387,5,0
{ 0 4 7 8 9 },-0.0216277,{ 0 1 4 7 8 },-0.0285387,5,1
{ 0 4 6 8 9 },-0.016636,{ 0 1 4 7 8 },-0.0285387,5,2
{ 7 8 },-0.0250616,{ 7 8 },-0.0250616,2,0
{ 3 8 },-0.0161735,{ 7 8 },-0.0250616,2,1
{ 0 4 7 8 },-0.0255352,{ 0 4 7 8 },-0.0255352,4,0
{ 0 1 7 8 },-0.0249667,{ 0 4 7 8 },-0.0255352,4,1
{ 0 1 8 9 },-0.0190089,{ 0 4 7 8 },-0.0255352,4,2
{ 2 3 5 6 },-0.0254904,{ 2 3 5 6 },-0.0254904,4,0
```
#### n10-doa.csv
A csv file of the minimal cost splits, their size & score, and the size of their domain of attraction.
```
Split,SplitSize,FinalScore,DoASize
1100100111,4,-0.0254904,110
1100100110,5,-0.0285387,105
1000100110,4,-0.0255352,75
1111011001,3,-0.025536,51
1111100110,3,-0.023584,35
1100111110,3,-0.0237436,22
1111111001,2,-0.0250616,17
1111100111,2,-0.0245457,13
1100111111,2,-0.0236405,7
1100000000,2,-0.022339,2
```

#### n10-graph.csv
A dot- (graphviz) format file producing this image:
![n10-graph-clust-anon.png](https://github.com/mcharleston/Flatbush/blob/master/n10-graph-clust-anon.png)

Alternatively with the defaults left as they are:
![n10-graph-defaults.png](https://github.com/mcharleston/Flatbush/blob/master/n10-graph-labelled.png)

### NEXUS format input files

Note that FlatBush can also parse fairly basic NEXUS-style input files.

The advantage of using a NEXUS(-like) input format is that users can also input sets of splits to analyse if they are of particular interest.
The input format is case IN-sensitive.

The fairly basic NEXUS-style format is as follows:

#### NEXUS file
```
	 <nexus> ::= "#NEXUS" { <nexusblock> }
	 <nexusblock> ::= "begin" <blockname> { <setting> } ( "end" | "endblock" ) ";"
	 <blockname> ::= <string>
	 <setting> ::= <id> [ "=" ] <val> [","]
```

#### Data Block
```
	 <datablock> ::= "begin" "flatbush" ';' { <component> } <endblock>
	 <component> ::= ( <dim> | <format> | <matrix> )
	 <dim>       ::= "dimensions" { ( <ntax> | <nchar> ) } ';'
	 <format>    ::= "format" { ( <datatype> | <missing> | <gap> ) } ';'
	 <matrix>    ::= "matrix" { ( <id> <sequence> ) } ;
	 <endblock>  ::= ( "end" | "endblock" ) ';'
```

#### Newick Tree format
```
	 Newick format grammar, adapted from http://evolution.genetics.washington.edu/phylip/newick_doc.html
	 
	 <tree>            ::= <descendant_list> [ <rootlabel> ] [ ':' <branchlength> ] ';'
	 <descendant_list> ::= '(' <subtree> { ',' <subtree> } ')'
	 <subtree>         ::= (
	                          <descendant_list> [ <internal_label> ] [ ':' <branchlength> ]
	                       |
	                          <leaf_label> [ ':' <branchlength> ]
	                       )
	 <rootlabel>       ::= <label>
	 <internal_label>  ::= <label>
	 <leaf_label>      ::= <label>
	 <label>           ::= <string>
	 <branchlength>    ::= <number>
	 
	 Further, the Newick format can be preceded by an optional tree name:
	 <treename>        ::= "tree" [ '=' ] <string>
```

#### Taxa Block
```
	<taxablock> ::= "taxa" ';' { <component> } ( "endblock" | "end" ) ';'
	<component> ::= <dim> | <taxlabels>
	<dim>       ::= "dimensions" { ( <ntax> | <nchar> ) } ';'
	<taxlabels> ::= "taxlabels" { <seqID> } ';'
```

#### FlatBush Block
```
  <flatbushblock> ::= "begin" "flatbush" ';' { <flatbushcomponent> } <endblock>
	<flatbushcomponent> ::= ( <tree> | <splits> | <outputformat> | <patternfrequencies> )
	<tree>         ::= <newicktree>
	<splits>       ::= "splits" ( <binarysplits> | <taxonsplits> ) ';'
	<binarysplits> ::= { <binarystring>  [ ',' ] }
	<taxonsplits>  ::= ( "taxonset" | "taxa" ) { <taxonset> [ ',' ] }
	<taxonset>     ::= '{' { <taxonname> [ ',' ] } '}'
	<outputformat> ::= ( "outputformat" | "outfmt" | "-of" ) (
	                          <fmt_binarycsv> |
	                          <fmt_binarytabbed> |
	                          <fmt_sets> ) ';'
	<fmt_binarycsv> ::= ( "0" | "binary_csv" )
	<fmt_binarytabbed> ::= ( "1" | "binary_tabbed" )
	<fmt_fmt_sets> ::= ( "2" | "sets" | "long" )
	<patternfrequencies> ::= ( "patternfrequencies" | "patternfreqs" | "pf" |
	                          "patternweights" | "pw" ) { ( <pattern> <number> ) [ ',' ] } ';'
	<endblock>     ::= ( "end" | "endblock" ) ';'
```
