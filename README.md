# Wild-pBWT

A PBWT-based algorithm for identifying all Maximal Perfect Haplotype Blocks with Wildcards (MPHBw).

## Build the tool(s)

### Prerequisites

A C++11 compiler and `make`. To read VCF/BCF files, [htslib](https://github.com/samtools/htslib) too (for example `sudo apt-get install libhts-dev` on Ubuntu, `brew install htslib` on macOS, or `conda install -c bioconda htslib`). It is found with `pkg-config`; without it everything else works and `make HTSLIB=0` leaves it out. If `pkg-config` does not know where it is, run `make HTSLIB=1 CPPFLAGS=-I<include dir> LDLIBS="-L<lib dir> -lhts"`.

### Build

Run
```
make
```
inside the main folder, to build `wild-pbwt`, `gen`, and `err`. You could also selectively build a tool by running
```
make <tool-name>
```
For example run `make wild-pbwt` to only build `wild-pbwt` tool.

## Usage

---

### `wild-pbwt`

The `wild-pbwt` binary under the `bin` subfolder computes all the maximal haplotype blocks with wildcards from a given input, with the extended pBWT

```sh
./bin/wild-pbwt -f <filename> -a <t-alleles> [-c|-o y|-r y] [-v y] [-b <block_size>] [-g <buffer_size>] [-w y] [-t y]
```
where:
- `<filename>` is the input file containing the haplotype panel in ASCII format where each line represents a single haplotype and each column is a variation site. Wildcards are represented with character `*`. All the lines must have the same number of characters (no spaces, no Windows `\r`), and a character that is neither `*` nor a digit below `-a` stops the run, saying in which haplotype and column it is. A file ending in `.vcf`, `.vcf.gz`, `.vcf.bgz` or `.bcf`, or `-` for a VCF/BCF on the standard input, is read as a VCF/BCF (see below).
- `<t-alleles>` is the alphabet size (not counting `*`), i.e., the maximum number of alleles in a single site. For a VCF/BCF it defaults to 2.
- `<block_size>` is the minimum block size required to count a maximal block. It defaults to 2. It requires the `-o` flag.
- `<buffer_size>` makes `wild-pbwt` read the input with one buffered file stream per haplotype (each with a buffer of `buffer_size` bytes, e.g. `512`) instead of mapping the whole file in memory, which is the default. Use it for panels that do not fit in memory. It needs one open file per haplotype.
- `-t y` says that the file is the matrix transposed: each line is a column (a site) with one character per haplotype, so that it is read in a single pass (see below).

A block is a set `K` of at least two haplotypes and a range of columns `i..j` such that the haplotypes of `K` do not disagree on any column of the range, the block cannot be extended to the left or to the right, and no other haplotype can be added to `K`. What "do not disagree" and "cannot be extended" mean for columns where all the haplotypes of `K` have a wildcard depends on the mode:

- By default `wild-pbwt` reports the maximal perfect haplotype blocks with wildcards of [Williams & Mumey](https://doi.org/10.1016/j.isci.2020.101149): on every column of a block the haplotypes of `K` hold exactly one allele besides wildcards, so a column with only wildcards can never be inside a block. A block cannot be extended over a column where its haplotypes hold two alleles, or only wildcards.
- With `-w y`, a column where all the haplotypes of `K` have a wildcard does not disagree with anything: a block goes on through it, and only a column with two alleles stops it. Each block is reported once.

In both modes one bit per cell of the panel is kept in memory (`M`x`N`/8 bytes) to look up wildcards.

#### VCF/BCF input

A VCF gives the panel one site at a time, which is what the pBWT needs, so it is read in a single pass with no conversion and no matrix on disk (`-g` has no meaning for it). Each record is a column, and the genotypes (`FORMAT/GT`) give the haplotypes: haplotype `h` of sample `s` is the row `ploidy * s + h`, and a missing allele (`.`) is a wildcard. The alleles are the numbers of the genotypes (`0` is `REF`, `1` the first `ALT`, and so on), so a site with more alleles than `-a` stops the run and says which `-a` to use. The genotypes should be phased (`|`): unphased ones are taken in the order they are written, with a warning. Blocks are reported with the row and the column indices, as for a matrix. Samples and positions are in the same order as in
```sh
bcftools query -l file.bcf                        # sample s is line s + 1
bcftools query -f '%CHROM\t%POS\n' file.bcf       # column c is line c + 1
```
All records are read, so select the sites before, for example biallelic SNPs of one region:
```sh
bcftools view -r 20:1-5000000 -v snps -m2 -M2 file.bcf | ./bin/wild-pbwt -f - -o y
```

#### Transposed matrix input

With `-t y` a matrix file is read by columns instead of by haplotypes: line `c` has the alleles of site `c`, one character per haplotype (`0`, `1`, ... and `*`, no separators), so the first line says how many haplotypes there are and all the lines must have as many characters. Like a VCF it is read in a single pass, with one open file and memory that does not depend on the number of sites, and `-f -` reads the standard input. Empty lines are skipped and a `\r` at the end of a line is ignored. A character that is not `*` or an allele below `-a` stops the run (`-a` is needed, as for any matrix). For example, the matrix
```
0101
0110
```
is the file `00`, `11`, `01`, `10` when transposed.

If ran with `-o` flag, blocks will be output to standard output. Blocks colud be saved to external file adding `> output_file.txt` at the end of the `wild-pbwt` command. A block is a line `[rows], first column, last column`, with rows and columns counted from 0, for example `[0,1,3,4,], 0, 2` for haplotypes 0, 1, 3 and 4 on columns 0 to 2. The rows are in no particular order.

If ran with `-r y` instead, the blocks are output in the same way but with the rows sorted and the runs of consecutive rows written as ranges, and without the final comma: `[0-1,3-4], 0, 2`. When many of the haplotypes of a block are neighbours in the panel, as in a panel with many similar haplotypes, this is much shorter (and quicker to write: on a panel of 500,000 haplotypes with 377,447 blocks the output went from 5.4 GB to 2.3 GB, and the run from 35 to 6 seconds).

If ran with `-c` flag, blocks will not be output, only counted.<br>
To specify a block_size, run with `-o` (`-c` not permitted).

If ran with `-v` flag, the extended pBWT execution will be output to the standard error.

---

### `gen`

The `gen` binary under the `bin` subfolder allows the user to generate a plain random matrix M (haplotypes) x N (SNPs). To add missing data to it, run `err` on the result.
```sh
./bin/gen <save_directory> <t-alleles> <haplotypes_count> <SNPs_count>
```

Example:
```sh
mkdir data
./bin/gen data 2 1000 52000
```
Will generate a matrix under the newly created `data` folder, `bi-allelic` with `1000` haplotypes and `52000` SNPs.

---

### `err`

The `err` binary under the `bin` subfolder allows the user to generate a matrix with missing data, from an input matrix. The rate of missing data insertion for each locus is specified by `wild_rate`\%.

```sh
./bin/err <input_matrix> <wild_rate> <path_to_output_matrix>
```
Example:
```sh
./bin/err data/input_matrix 3 data/output_matrix 
```
Will generate `output_matrix` under `data` folder as a copy of `input_matrix` with 3\% missing data.


## Testing

`test/check_blocks.py` compares the blocks found by `wild-pbwt` with a brute force on random small panels, and with a slower but independent search on a single panel of any size. Both follow the definitions above, and the exit code is not 0 if a block is missing, reported twice, or not a block:

```sh
python3 test/check_blocks.py bin/wild-pbwt 2000 1 wm                         # default mode
python3 test/check_blocks.py bin/wild-pbwt 2000 1 wild                       # -w y
python3 test/check_blocks.py bin/wild-pbwt --panel data/<panel> -a 2 wm      # one panel
```
