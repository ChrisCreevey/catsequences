# catsequences

### A tool for concatenating multiple fasta alignments for supermatrix phylogenetic analyses

Please reference this DOI: [![DOI](https://zenodo.org/badge/79135872.svg)](https://zenodo.org/badge/latestdoi/79135872)


Usage: `catsequences FILENAME`

Where:
  - `FILENAME` contains a list of files with sequences to be concatenated.
  - The concatenated alignments are outputted to a file called `allseqs.fas`.
  - Partition information will be written to the file `allseqs.partitions.txt`

The sequence names are read up to the first `.` `|` or ` ` (space) or the end of the name.

**NOTE**: You can change the delimiters that are used to split the name by editing the code - [see below for details](https://github.com/ChrisCreevey/catsequences/edit/master/README.md#i-need-a-different-delimiter-character-in-the-sequence-names-than-the-default---or---space)


----------------------------------

# To install:

Download the file `catsequences.c` and or if you have git installed use the command:
```
git clone https://github.com/ChrisCreevey/catsequences.git
```

To build catsequences type:
```
cc catsequences.c -o catsequences -lm
```

Copy the executable `catsequences` to somewhere on your path like `~/bin` to have access to it from anywhere.

--------------------------------

# To use:

Usage: 

```
catsequences FILENAME
```

  where `FILENAME` is a file containing a list of files with sequences to be concatenated.

  - The concatenated alignments are outputted to a file called `allseqs.fas`
  
  - Partition information will be written to the file `allseqs.partitions.txt`

  - The sequence names are read up to the first `.`, `|` or ` ` (space) (or the end of the name)

The easiest way to create the list input file is to use a command like:
```
ls *.fasta > list
```
assuming that all your input file have the extension ".fasta"

The software will put in '?' for the entire length of any gene for which the species is missing.


# Trouble Shooting

## If the total concatenated seuquences length is larger than defined in the code you recieve an ERROR saying:

```
ERROR: The total concatenated sequence length is greater than the maximum length defined (1000000)!
    To fix this increase the value of MAXLEN on line 15 of 'catsequences.c` and recompile the tool.
```

This line you want to change is:

```
#define MAXLEN 1000000
```

Although this means that the total length of all the seqeunces you are trying to concatentate is greater than 1 million bp. 
Please make sure that this is as expected before changing the code.

## I need a different delimiter character in the sequence names than the default "|", "." or " " (space)

Add (or remove) characters between the `" "` as needed on line 16 of `catsequences.c`:

```
#define NAMEDELIMITER " |." 
```

Then recompile (build) the tool.


