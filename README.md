# PanDelos plus plus
PanDelos plus plus: a dictionary-based method for pan-genome content discovery

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [](#lang-en)

----

## Briefly description

**PanDelos plus plus** is used for the discovery of the pangenomic content in the population of bacteria with a type of approach without alignment and based on the analysis of the multiplicity of k-mer.
It is a natively scalable methodology, whose algorithms are executed in parallel with OpenMP.
----

## Software architecture
The PanDelos plus plus software is organized in 2 python modules and a C++ software that are piped together by a bash script.
The bash script, `pandelos.sh`, provides the acces point the the PanDelos plus plus pipeline.

```
bash pandelos.sh <input.faa> <output_graph.net> <sequences_type> <log_file.txt> <clus_file.clus> <coco_file.txt>
```

<hr />

### Input format

The complete set of (gene) sequences `<input.faa>`, belonging to any of the studied genomes, must be provided as a text file.


For each sequence, two lines are reported in the file. An identification line that is composed of three parts separated by a tabulaiton character. The parts represent the genome identifiers, the gene identifier and the gene product.

After the identification line, the complete gene sequence in FASTA amino acid format is reported in a single line. No black lines are admitted between the indetification line and the sequence line, neighter between genes.

A valid file is given by the following example listing 4 genes from 2 genomes:
```
NC_000913	b0001@NC_000913:1	thr operon leader peptide
MKRISTTITTTITITTGNGAG
NC_000913	b0024@NC_000913:1	uncharacterized protein
MCRHSLRSDGAGFYQLAGCEYSFSAIKIAAGGQFLPVICAMAMKSHFFLISVLNRRLTLTAVQGILGRFSLF
NC_002655	Z_RS03160@NC_002655:1	hok/gef family protein
MLTKYALVAVIVLCLTVPGFTLLVGDSLCEFTVKERNIEFRAVLAYEPKK
NC_002655	Z_RS03165@NC_002655:1	protein HokE
MLTKYALVAVIVLCLTVLGFTLLVGDSLCEFTVKERNIEFKAVLAYEPKK
```
> IMPORTANT: make sure the gene identifiers are unique within the input file. Commonly used file formats used to share genome annotaitons do not require that different locus tags of the same gene must be unique.

We suggest to use the following format to build unique gene identifiers:
```
gene_identifier@genome_identifier:unique_integer
```
The fields `gene_identifier` and `genome_identifier` are the same reported in the input file, while the `unique_integer` is used to disitrnghuish multiple copies of the same gene (same gene identifier) wihtin the same genome. The integer starts from 1 and it is incremented according to the order gene are written in the input file.


The examples provided in the `examples` folder generate 4 different dataset files, having the `.faa` extension, which can be consulted.

<hr />

### Output graph
Running the module in C ++ produces a graph named `[output_graph].net` containing the gene families, within which inconsistent families may be found

### Sequences type
This software allows the discovery of pangenomic content for nucleotide (1) or amino acid (0) sequences

### Log file
The execution requires a log file which will contain some main information about the execution

### Clus file
The execution of PanDelos produces an output file named `[clus_file].clus` which reports the gene families retrieved by the software.
Each row of the output file represented a specific gene family retrieved by PanDelos.

<hr />

## Installation
### Requirements
Before running ***PanDelos plus plus***, please verify that the following software is installed on your system
* bash
* zip
* g++ version 6 or higher
* C++17 or higher
* OpenMP 5.2 or higher
* GLIBCXX 3.4.2 or higher
* python version 3.7 or higher
* biopython
* bcbio-gff ( https://github.com/chapmanb/bcbb/tree/master/gff )

### Setup
- Download the software from here or clone the github repository
- Enter the `PanDelos-plusplus` directory and type
```
bash setup/setup.sh
```
to compile the C++ source code of ***PanDelos-plus-plus*** and create the folders that the software needs.

The script allows to download from a repository of the files concerning genomes generated with PANPROVA of Mycoplasma genitalium and Escherichia coli.

These files, using utility scripts that will be described later, will be transformed into .faa files

Furthermore, it is also possible to have the same genome files in the form of nucleotides (for a total of 4 datasets available)

### Compiling the C++ source code
The script setup/compile.sh is available to manually compile the program in C ++ or you can give the command:

```
g++ -w -std=c++17 -fopenmp -O3 src/cpp/main.cpp -o bin/pandelos_plus_plus.out

```
The software needs the C ++ 17 standard and the -fopenmp directive which enables the compiler to manage the pragma directives of the library that deals with parallelizing the algorithms in the software

<hr />

## Running the examples
Examples are available in the folder examples/ o to test PanDelos plus plus.

In general, each example script allows you to start multiple tests, separate and sequential, involving a different number of genomes, even in the form of intervals.

There are 4 examples available:
1) 16 genomes of Mycoplasma genitalium with nucleotide sequences
2) 16 genomes of Mycoplasma genitalium with amino acid sequences
3) 16 genomes of Escherichia coli with nucleotide sequences
4) 16 genomes of Escherichia coli with amino acid sequences

----

## Utilities

* `script/panprova2gbk.sh.py`: a Python script for converting a GFF file into the GBK file
* `script/panprova2nucleotides.sh script/panprova2nucleotides.sh`: a Python script for converting a FNA and GFF files into a FASTA file composed of nucleotide sequences.
* `script/gbk2faa.sh gbk2faa.py`: a python script to convert one or more GBK files into a single FASTA file
* `src/python/quality.py`: calculates statistics about th eextracted pan-genome content and print them.

----

## License
PanDelos plus plus is distributed under the MIT license. This means that it is free for both academic and commercial use. Note however that some third party components in PanDelos plus plus require that you reference certain works in scientific publications.
You are free to link or use PanDelos plus plus inside source code of your own program. If do so, please reference (cite) PanDelos plus plus and this website. We appreciate bug fixes and would be happy to collaborate for improvements.

<hr />

<!--
## Detailed description

The following picture gives a detailed description of the PanDelos plus plus workflow.

<p align="center">
<img src="https://github.com/InfOmics/PANPROVA/blob/main/workflow.svg?raw=true" alt="workflow" width="400"/>
 </p>

The workflow is composed of a set of internal tools, Python scripts and C++ executables, plus some external Python scripts that can be used for file format conversions.

Sections with a yellow background are those internal tools that are in charge of the `PANPROVA.sh` script.

<br/>

The internal tools are:
* `create_hgt_pool`: a C++ executable for creating an HGT pool from a set of PEG files. It also takes are input the root genome in roder to discard genes that a re similar to the genetic sequences within the root genome.
* `generate_tree.py`: a Python script for randomly generating a phylogenomic tree of the wanted population.
* `tree2phyloxml.p`: a tool for converting a PANPROVA tree into a PhyloXML file and for generating an image showing it.
* `evolve`: a C++ executable that implements the evolution procedure.
* `get_pan_distrs.py`: a Python script for retrieving pangenomic information from the generated population and for creating the corresponding output.
* `pegs2gxx.py`: a Python script for converting the generated genomes into the GBK and GFF+FASTA formats.

----

### Extraction of HGT pool

The pool of HGT genes to be used during the evolution simulation is extracted from a set of input genomes (in PEG format) and by taking into account genes that are already present in the root genome (still in PEG format) for excluding genes similar to them from the HGT pool.
The following picture illustrates the main steps of the extraction procedure.

<p align="center">
<img src="https://github.com/InfOmics/PANPROVA/blob/main/createhgt.svg?raw=true" alt="create hgt" width="200"/>
 </p>

From the given input genomes, a set of genes that are not similar to the genes present in the root genome is initially extracted. Then a nonredundant pool of genes is created by discarding genes that are similar to other genes in the initial set.
The similarity among nucleotide genetic sequences is computed by taking into account the similarity between their k-mer content [1]. In particular, a Jaccard similarity between k-mer multisets of two genetic sequences is computed. Genes with a similarity greater than 0.3 with root genes are discarded. Successively, we set an arbitrary order of the surviving genes. Then, each gene is compared with genes that come after it in the ordering. If the similarity is greater than 0.5, then the latter gene is marked to be discarded. At the end of the scanning, all the genes that were marked are removed from the HGT pool.

### Evolution procedure

The workflow of the evolution procedure, together with examples (in yellow boxes) of intermediate data, is shown in the following figure.

<p align="center">
<img src="https://github.com/InfOmics/PANPROVA/blob/main/evolve.svg?raw=true" alt="evolve" width="500"/>
 </p>

The workflow refers to the case in which the generation of the random phylogenomic tree is integrated into the process.
<br/>

At each step, a genome from the current population is chosen to be the parent of the next genome to be created. Thus, the parent genome is cloned and an initial version of the child genome is produced (see example 1 of the figure).
<br/>

Then, according to a given probability, each vertically transmitted gene is selected to be altered or not. If yes, its loci are variated according to a given variation percentage. Possible variations are substitution, insertion or deletion.
The tool gives the possibility to specify user-defined substitution probabilities for nucleotides by providing a file containing them. By default, every nucleotide can be substituted by any other nucleotide with equal probability.
Any modification is applied such that it does not produce or modify any star or stop codon of genes that overlap the gene that is currently modified. Overlapping genes may reside on both strands.
Because valid genetic sequences must be provided, substitution regards one nucleotide at a time, while insertion and deletion regard 3 nucleotides at a time, such that the length of the resulting sequence is still a multiple of 3.
<br/>
Ts/Tv ratio and synonym/non-synonym mutation ratio are intended to be the effects of the alterations that are performed on genetic sequences, thus they can not be specified as input parameters. We are aware that more complex models of sequencing alteration are available at the state of the art. However, the main aim of  ***PANPROVA*** is to simulate pangenomic effects, mainly due to the acquisition and deletion of genes. An extension of the software by us or the research community may include more accurate models.
<br/>

Subsequently, variated vertically transmitted genes are selected to be duplicated within the new genome according to a given probability.
<br/>

Duplication, insertion of HGT genes and transposition of genes is made such that a random locus of the genome is chosen. the locus must not be covered by any other gene. Thus, the genetic sequence of the gene, together with start and stop codons, is inserted at the selected locus. See examples 2 and 4 of the figure.
The resultant gene set is modified by a given percentage. If the set is composed of n genes and 2% of the set has to be variated, then (n/100)x2 variation operations are performed. such operation can be a horizontal gene acquisition of a gene removal. If the probability that an operation is an acquisition is p, then the probability that the operation is a removal is 1-p.
<br/>

In the case of gene removal, a gene is randomly chosen to be removed. All the nucleotides that belong to the selected gene are removed from the genome if they do not overlap other genes. See example 3 of the figure.
<br/>

In case of gene acquisition, if the HGT pool is not empty, a genetic sequence is randomly chosen from the pool, inserted in the genome and removed from the pool. See example 4 of the figure. If the HGT pool is empty, a purely random nucleotide sequence is generated and inserted within the genome.
<br/>

Subsequently, the resultant set of genes is randomly picked for transposition according to a given probability.
<br/>

Lastly, the new genome is added to the population and the process is repeated until the desired number of genomes is produced. Every time a new genome is produced, its parenting relationships are recorded. In particular, the information regarding the genome from which it has been cloned is stored. In addition, for each gene in the new genome, the information regarding the parent gene is stored. for vertically transmitted genes, such information reports the identifiers of the gene present in the parent genome. For duplicated genes, such information reports the identification of the paralog gene from which the gene has been duplicated. For horizontally transmitted genes, such information is null. See example 5 of the Figure.

