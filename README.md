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
usage: flatbush [OPTIONS]
	-f <int>: output format for errors and splits.
		0:	binary representation of each split, comma-separated-values (CSV);
		1:	binary representation of each split, tab-delimited;
		2:	splits as sets, comma-separated (CSV).
		Default: 0 (binary, csv)
	-h:	 stream this usage and quit
	-i <filename>: input file name in Fasta or Nexus format.
		No default value --- without this, nothing much happens.
	-l: do a landscape analysis.
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
	--sdrs: do a single Steepest Descent Random Start (useful for debugging but not much else)
	--seed <int>: set the random number seed, for repeatability
		Default: seed set from the system clock
	--show-all-details: show all the split errors:
		this will potentially take a much longer time and fill up your terminal with junk you didn't really want.
	--silent flag to turn on Silent mode, useful for bash scripting
		Default: off
		--splits-file <filename>: load a set of splits in as "true" splits
		format: a list of splits as binary strings with a floating-point number for each, such as branch length.
		e.g. "00111, 0.5" to define split {0,1} | {2,3,4} with value 0.5
	--taxon-names: if possible, print splits with taxon names;
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

Note that FlatBush can also parse fairly basic NEXUS-style input files.

Another example:
```
> ./FlatBush -i noclock-GTR-n20-b0.5-l1000-1-2.fst -m clust -n noclock-GTR-n20-b0.5-l1000-1-2.fst.clust --splitsfile noclock-GTR-n20-b0.5-l1000-1.bsplits -l
```
also loads a *splits* file, containing a set of splits of interest.
[more to come here]


### NEXUS format input files

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
