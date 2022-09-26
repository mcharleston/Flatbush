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
	-h or no arguments at all: print this information.
	-i <filename>: input file name in Fasta or NEXUS format.
		No default value --- without this, nothing much happens.
	-l: do a landscape analysis.
	-o <filename>: output file name for results.
		Default value is "flatbush.out".
	-f <int>: output format for errors and splits.
		0:	binary representation of each split, comma-separated-values (CSV);
		1:	binary representation of each split, tab-delimited;
		2:	splits as sets, comma-separated (CSV).
	[-n|-name] <string>: a project name, from which output files will be derived automatically.
	-pert <int>: which perturbation types to be used for navigating the split space.
		1: single taxon move only (move a single taxon from one side of the split to the other);
		2: taxon swap only (swap taxa between the two sides of the split, maintaining the size of each part);
		3: allow both perturbations as for (1) and (2) above.
	-hsv <float> <float> <float>: input hue, saturation and value for the colouration of output graphs.
		(NB: saturation is in fact ignored: still trying to find something that makes sense to vary.)
	-s flag to turn on Silent mode for scripting:
		Default value is FALSE.
	-ss <int>: set the sample size of starting splits for landscape analysis.		Default value: 1000.
	--do-all-splits: analyse all splits (not usable if the number of taxa is more than 16).
	--sdrs: do a single steepest descent from a random starting split and then exit.
		Default value: FALSE.
	--seed <int>: set the random number seed.
		Default value: set by system clock.
	--show-all-details: output a potentially massive quantity of extra details that you probably don't want.
		Default value: FALSE.
	-v flag to turn on Verbose mode:
		Default value is FALSE.
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
