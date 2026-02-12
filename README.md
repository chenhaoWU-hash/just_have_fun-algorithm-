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
### 1.Creating an index
Generate SPSS and FM-index using FASTA files as input:  
```
python3 SequencesToIndexedSPSS.py -i <input_file> -k <kmer_size> -t <threshold> -o <output_file>  
```
-**Available options**
```
-i : the input FASTA file.  
-k : k-mer size (positive integer)).  
-t : filtering threshold (positive integer)).  
-o : output .dump file (FM-index).
```
-**Exemple** :   
```
python3 SequencesToIndexedSPSS.py -i reads_30x.fasta -t 2 -k 31 -o my_fm_index
```
### 2.Executable query ### 
Comparison of query sequence similarity using the generated FM-index file:  
```
python3 query_indexed_spss.py -q <query_file> -i <index_file> -k <kmer_size> -o <output_file>  
```
## Structure of the project ##
```

|-- README.md             # Document de description du projet
|-- Code/
|   |-- query_indexed_spss.py       # Module d'interrogation
|   |-- fmi.py                      # FM-index réalisé
|   |-- SequencesToIndexedSPSS.py # Module de construction d'index
|   `-- timer.py                    # mesure du temps
|   |-- spss.py                     # Module de construction de spss
|   |-- tools_karkkainen_sanders.pu  # calcule de SA pour FM-index
```

