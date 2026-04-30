# catsequences

### A tool for concatenating multiple FASTA alignments for supermatrix phylogenetic analyses

Please reference this DOI: [![DOI](https://zenodo.org/badge/79135872.svg)](https://zenodo.org/badge/latestdoi/79135872)

---

## Installation

Download `catsequences.c` or clone the repository:

```bash
git clone https://github.com/ChrisCreevey/catsequences.git
```

Build:

```bash
cc catsequences.c -o catsequences -lm
```

Copy the executable to somewhere on your `PATH` (e.g. `~/bin`).

---

## Basic usage

```
catsequences [OPTIONS] LISTFILE
catsequences [OPTIONS] FILE1.fasta FILE2.fasta ...
```

**List-file mode** (original behaviour):

```bash
ls *.fasta > list
catsequences list
```

**Direct FASTA mode** (pass files on the command line):

```bash
catsequences gene1.fasta gene2.fasta gene3.fasta
```

Glob wildcards are also accepted inside a list file:

```
genes/*.fasta
```

---

## Options

| Flag | Description |
|------|-------------|
| `-d DELIMS` | Name-delimiter characters (default: `" \|."`) — truncate the species name at the first occurrence of any character in this string |
| `-o PREFIX` | Output filename prefix (default: `allseqs`) |
| `-m CHAR` | Missing-data character (default: `?`; common alternatives: `-`, `N`) |
| `-q` | Quiet — suppress all progress messages |
| `-v` | Verbose — print per-file scan/read progress |
| `-f FORMAT` | Partition file format: `generic` (default), `raxml`, `iqtree`, `nexus` |
| `-F FORMAT` | Alignment output format: `fasta` (default), `nexus`, `phylip` |
| `-O` | Write an occupancy matrix CSV |
| `--min-taxa N` | Skip any gene file where fewer than `N` × 100 % of species are present (e.g. `0.5` = require at least 50 % occupancy) |
| `--on-duplicate MODE` | What to do when a species name appears more than once in a file: `warn` (default — keep first, print a message), `error` (abort), `keep-last` (overwrite with the last copy) |
| `--include FILE` | Only output species whose names appear in `FILE` (one name per line, without `>`) |
| `--exclude FILE` | Exclude species whose names appear in `FILE` (one name per line, without `>`) |

---

## Outputs

| File | Description |
|------|-------------|
| `<PREFIX>.fas` / `.nex` / `.phy` | Concatenated alignment |
| `<PREFIX>.partitions.txt` | Partition coordinates |
| `<PREFIX>.occupancy.txt` | Per-species and overall occupancy summary |
| `<PREFIX>.occupancy.csv` | Occupancy matrix — present (1) / absent (0) per gene per species (written only with `-O`) |

Species that are absent from a gene file are padded with the missing-data character for that gene's full length.

---

## Partition format examples

### generic (default)
```
gene1.fasta=1-500;
gene2.fasta=501-1200;
```

### raxml (`-f raxml`)
```
DNA, gene1 = 1-500
DNA, gene2 = 501-1200
```

### iqtree (`-f iqtree`)
```
gene1 = 1-500
gene2 = 501-1200
```

### nexus (`-f nexus`)
```
#NEXUS
begin sets;
charset gene1 = 1-500;
charset gene2 = 501-1200;
end;
```

---

## Examples

Run with a list file, writing RAxML partitions and a NEXUS alignment:

```bash
catsequences -f raxml -F nexus -o mymatrix list
```

Require at least 50 % of species per gene, write occupancy CSV:

```bash
catsequences --min-taxa 0.5 -O -o filtered list
```

Include only a subset of species:

```bash
echo -e "Homo_sapiens\nPan_troglodytes\nGorilla_gorilla" > primates.txt
catsequences --include primates.txt -o primates list
```

Change the delimiter so only `|` splits names (not `.` or space):

```bash
catsequences -d '|' list
```

---

## Testing

Build and run the regression test suite:

```bash
cc catsequences.c -o catsequences -lm
bash tests/run_tests.sh
```

---

## Notes

- Species name lengths and the number of species and genes are all determined
  dynamically at run time — there are no hard-coded limits requiring recompilation.
- The tool reads each input file twice: once to discover species and measure
  sequence lengths, and once to read sequence data.  This keeps peak memory
  proportional to `(number of species) × (total concatenated length)`.
