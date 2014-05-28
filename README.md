Preprocessor
============

A Ruby gem for preprocessing mRNA reads from Illumina sequencing.

# Dependencies

 - Ruby
 - Java
 - Khmer
 - SPAdes

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