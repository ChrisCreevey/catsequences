# catsequences

### A tool for concatenating multiple fasta alignments for supermatrix phylogenetic analyses


Usage: `catsequences filename`

Where filename contains a list of files with sequences to be concatenated
The concatenated alignments are outputted to a file called "allseqs.fas"
Partition information will be written to the file "allseqs.partitions.txt"

The sequence names are read up to the first "." or "|" (or the end of the name)


----------------------------------

# To install:

Download the file "catsequences.c" and or if you have git installed use the command:
```
git clone https://github.com/ChrisCreevey/catsequences.git
```
To build catsequences type:

```
cc catsequences.c -o catsequences -lm
```
Copy the executable "catsequences" to somewhere on your path like "~/bin" to have access to it from anywhere.

--------------------------------

# To use:

Usage: 

```
catsequences filename
```

  where filename contains a list of files with sequences to be concatenated

  >The concatenated alignments are outputted to a file called "allseqs.fas"
  
  >Partition information will be written to the file "allseqs.partitions.txt"

  >The sequence names are read up to the first "." or "|" (or the end of the name)

The easiest way to create the list input file is to use a command like:

```
ls *.fasta > list
```
assuming that all your input file have the extension ".fasta"

The software will put in '?' for the entire length of any gene for which the species is missing.

