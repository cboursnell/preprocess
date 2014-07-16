Preprocessor
============

A Ruby gem for preprocessing mRNA reads from Illumina sequencing.

Caveat: currently only works with paired reads.

# Input

The input is file containing one line per fastq file with 4 fields separated by commas. The 4 columns are:

 - Filename, preferably absoluate location not relative
 - Repetition. The sample number (Integer)
 - Type. This could be tissue name or cell type etc (String)
 - Pair. Either 1 or 2 for which file in the pair (Integer)

For example:

```
/home/chris/documents/rice/rice_BS_1_1.fq,1,BS,1
/home/chris/documents/rice/rice_BS_1_2.fq,1,BS,2
/home/chris/documents/rice/rice_BS_2_1.fq,2,BS,1
/home/chris/documents/rice/rice_BS_2_2.fq,2,BS,2
/home/chris/documents/rice/rice_M_1_1.fq,1,M,1
/home/chris/documents/rice/rice_M_1_2.fq,1,M,2
/home/chris/documents/rice/rice_M_2_1.fq,2,M,1
/home/chris/documents/rice/rice_M_2_2.fq,2,M,2
```

Currently only paired reads are supported.

# Usage

From the command line

```
preprocess --input <input> --output <output> --threads <t> --memory <m>
```

for example:

```
preprocess --input data --output ~/rice/output --threads 8 --memory 20 --verbose
```

# Dependencies

 - Ruby (at least version 2.0.0)
 - Java (at least -version 1.7.0_55)
 - Khmer (normalise-by-median.py)
 - SPAdes (automatically downloaded if not already installed)
 - trimmomatic (jar file included included)

# Installation of Khmer

First get the python development libraries and pip

```
sudo apt-get install python2.7-dev python-virtualenv python-pip gcc
```

Then

```
pip install khmer
```

# Installation of SPAdes

Download the compiled binaries

```
wget http://spades.bioinf.spbau.ru/release3.0.0/SPAdes-3.0.0-Linux.tar.gz
tar -xzf SPAdes-3.0.0-Linux.tar.gz
```