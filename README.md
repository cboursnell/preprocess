Preprocessor
============

A Ruby gem for preprocessing mRNA reads from Illumina sequencing.

Caveat: currently only works with paired reads.

# Input

The input is file containing one line per fastq file with 4 fields separated by commas. The 5 columns are:

 - Experiment name
 - Filename, preferably absoluate location not relative
 - Repetition. The sample number (Integer)
 - Type. This could be tissue name or cell type etc (String)
 - Pair. Either 1 or 2 for which file in the pair (Integer)

For example:

```
rice,/home/chris/documents/rice/rice_BS_1_1.fq,1,BS,1
rice,/home/chris/documents/rice/rice_BS_1_2.fq,1,BS,2
rice,/home/chris/documents/rice/rice_BS_2_1.fq,2,BS,1
rice,/home/chris/documents/rice/rice_BS_2_2.fq,2,BS,2
rice,/home/chris/documents/rice/rice_M_1_1.fq,1,M,1
rice,/home/chris/documents/rice/rice_M_1_2.fq,1,M,2
rice,/home/chris/documents/rice/rice_M_2_1.fq,2,M,1
rice,/home/chris/documents/rice/rice_M_2_2.fq,2,M,2
```

Currently only paired reads are supported.

# Usage

From the command line

```bash
preprocess --input <input> --output <output> --threads <t> --memory <m>
```

for example:

```bash
preprocess --input data --output ~/rice/output --threads 8 --memory 20 --verbose
```
or
```bash
preprocess --left reads_1.fq --right reads_2.fq --output ~/rice/output --threads 8 --memory 20 --verbose --trimmer none --correction hammer
```

# Dependencies

 - Ruby (at least version 2.0.0)
 - Java (at least version 1.7.0_55)
 - Khmer (normalise-by-median.py)
 - SPAdes (automatically downloaded if not already installed)
 - trimmomatic (jar file (automatically downloaded)
 - bbmap (automatically downloaded if not already installed)
 - skewer (automatically downloaded if not already installed)
 - facs (automatically downloaded if not already installed)


# Installation of Ruby

We recommend using the RVM to install Ruby. Installing with `sudo` can cause problems with certain gems.

Just run this to install the latest version of the RVM and Ruby

```bash
\curl -sSL https://get.rvm.io | bash -s stable --ruby
```

Then to install Preprocessor

```
git clone git@github.com:cboursnell/preprocess.git
cd preprocess
gem build *spec
gem install *gem
```


# Installation of Khmer

If you choose to use Khmer for digital normalisation you can install it using the python package manager pip. First get the python development libraries and pip

```bash
sudo apt-get install python2.7-dev python-virtualenv python-pip gcc
```

Then

```bash
pip install khmer
```
This might have to be run using `sudo`. This might not work at all... python package management is ... fraught. (https://twitter.com/gardaud/status/357638468572151808)

# Installation of SPAdes

BayesHammer is automatically installed when the program is run.

If you choose to use BayesHammer for read error correction it is included in the SPAdes assembler. The compiled binaries can be downloaded:

```bash
wget http://spades.bioinf.spbau.ru/release3.5.0/SPAdes-3.5.0-Linux.tar.gz
tar -xzf SPAdes-3.5.0-Linux.tar.gz
```

# Installation of BBMap

bbmap is automatically installed when the program is run.

# Installation of skewer

skewer is automatically installed when the program is run.
