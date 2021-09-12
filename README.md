# iPrimer
This repository holds the scripts used to design LAMP primers from an input fasta file.


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

* This package is supported for *Linux* operating systems.  The package has been tested on the following systems:
```
   Linux: CentOS Linux release 7.8.2003 (Core)
   Linux: Red Hat Enterprise Linux WS release 4 (Nahant Update 5)
   Linux: Ubuntu 18.04.4 LTS
```
* Perl 5 interpreter or higher is required.
   * [Installation instruction](https://learn.perl.org/installing/)
   * [Perl download](https://www.perl.org/get.html)
* Following tool packages placed in main folder are also required.
   * /NCBI/blast
   * /NCBI/dust

### Installation of iPrimer standalone program

* Place the iPrimer.tar.gz file anywhere in your Linux system and uncompress using the following command:
```
   tar -xzvf iPrimer.tar.gz   
```
* Copy your input FASTA files into the newly created iPrimer directory.
* Type 'perl iPrimer.pl' to run the program and view the help file.

### Command Line Parameters

* FASTA file submission, required (-f path/to/fasta).
   This option allows the user to submit one sequence in a FASTA file, using the following command:
```
   perl iPrimer.pl -f test.fasta
```
   This file should be provided in FASTA format.  In a FASTA file, a definition line that begins with begins with ‘>’ is required for each DNA sequence.  For example:
```
   >NC_045512.2:28274-29533 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
   ATGTCTGATAATGGACCCCAAAATCAGCGAAATGCACCCCGCATTACGTTTGGTGGACCCTCAGATTCAA
   CTGGCAGTAACCAGAATGGAGAACGCAGTGGGGCGCGATCAAAACAACGTCGGCCCCAAGGTTTACCCAA
   TAATACTGCGTCTTGGTTCACCGCTCTCACTCAACATGGCAAGGAAGACCTTAAATTCCCTCGAGGACAA
   GGCGTTCCAATTAACACCAATAGCAGTCCAGATGACCAAATTGGCTACTACCGAAGAGCTACCAGACGAA
   TTCGTGGTGGTGACGGTAAAATGAAAGATCTCAGTCCAAGATGGTATTTCTACTACCTAGGAACTGGGCC
   AGAAGCTGGACTTCCCTATGGTGCTAACAAAGACGGCATCATATGGGTTGCAACTGAGGGAGCCTTGAAT
   ACACCAAAAGATCACATTGGCACCCGCAATCCTGCTAACAATGCTGCAATCGTGCTACAACTTCCTCAAG  
```
   one fasta file containing a single sequence (1,260 nt) is provided for testing purpose.
* Apply relaxed condition, optional (-r y or -r yes).
   This option allows the user to design LAMP primers in relaxed condition.
```
   perl iPrimer.pl -f test.fasta -r y   
```
### Outputs

If the file is read in correctly, the following output files will be generated in iPrimer folder.

* Designed LAMP primers.  The results are made available in a tab-delimited text file.
```
   primerOutput.xls
```
* Log file contain iPrimer version, runtime, and input information.
```
   formatdb.log
```
## License & copyright

iPrimer is distributed under [Apache 2.0 with Commons Clause](LICENSE) license.