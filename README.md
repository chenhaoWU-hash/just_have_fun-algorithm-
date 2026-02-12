# just_have_fun-algorithm-
## Project description ##
We aim to represent and index a DNA sequencing file in order to quickly compare query sequences to this set of sequences. The comparison is performed without alignment and consists of counting the number of -mers shared between the query sequence and the sequences in the sequencing file.
## Functionnal ##
### 1. Generation of the SPSS index (Shared k-mers) ###
Creattion of an SPSS sequence containing all k-mers shared between the input sequences. Filtering based on a robustness threshold, ensuring greater accuracy of results. Supported input format: FASTA files.  
### 2.Construction of FM index ###
Using the Burrows-Wheeler Transform (BWT) algorithm and suffix arrays to construct an efficient FM-index. This index is designed to optimize pattern finding in large datasets.  
### 3. Searching and comparing sequences ###
Analysis of the similarity between a query sequence and the FM-index via the proportion of shared k-mers.
## Prerequisites ##
-Python 3.7+  
-package: pickle，argparse  
## How to use ## 
Generate SPSS and FM-index using FASTA files as input:  
```
python3 SequencesToIndexedSPSS.py -i <input_file> -k <kmer_size> -t <threshold> -o <output_file>  
```
## Available options ##
```
-i : the input FASTA file.  
-k : k-mer size (positive integer)).  
-t : filtering threshold (positive integer)).  
-o : output .dump file (FM-index).
```
-**Exemple** :   

