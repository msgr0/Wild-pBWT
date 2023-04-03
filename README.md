# Wild-pBWT

A PBWT-based algorithm for identifying all Maximal Perfect Haplotype Blocks with Wildcards (MPHBw).

## Build the tool(s)

### Prerequisites

SDSL by @SimonGog is required to build the project. More information available at [GitHub Repo](https://github.com/simongog/sdsl-lite).

You can either install it system-wide with the package manager (for example with `sudo apt-get install libsdsl-dev` on Ubuntu and similar), or using a conda environment (for example with `conda create --prefix ./local-env --channel conda-forge --override-channels compilers make sdsl-lite`, albeit [mamba on miniforge](https://github.com/conda-forge/miniforge) is suggested), or installing SDSL manually and modifying the `CXXFLAGS` in the `Makefile` (not suggested).

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
./bin/wild-pbwt -f <filename> -a <t-alleles> [-c|-o y] [-v y] [-b <block_size>] [-g <buffer_size>]
```
where:
- `<filename>` is the input file containing the haplotype panel in ASCII format where each line represents a single haplotype and each column is a variation site. Wildcards are represented with character `*`.
- `<t-alleles>` is the alphabet size (not counting `*`), i.e., the maximum number of alleles in a single site
- `<block_size>` is the minimum block size required to count a maximal block. It defaults to 2. It requires the `-o` flag.
- `<buffer_size>` is the buffer_size of the file stream buffer. No need to tweak this parameter in a normal use case. Adjusting this value according to the input column size could result in a performance improvement.

If ran with `-o` flag, blocks will be output to standard output. Blocks colud be saved to external file adding `> output_file.txt` at the end of the `wild-pbwt` command.

If ran with `-c` flag, blocks will not be output, only counted.<br>
To specify a block_size, run with `-o` (`-c` not permitted).

If ran with `-v` flag, the extended pBWT execution will be output to the standard error.

---

### `gen`

The `gen` binary under the `bin` subfolder allows the user to generate a matrix M (haplotypes) x N (SNPs) with a specified wildcard rate.
```sh
./bin/gen <save_directory> <t-alleles> <haplotypes_count> <SNPs_count> <wild_rate> 
```
If no `wild_rate` is given, a plain random matrix will be generated.

If a `wild_rate` is given, a plain random matrix will be generated alongside a version of that matrix with missing data.

Example:
```sh
mkdir data
./bin/gen data 3 1000 25000 5
```
Will generate 2 matrices under the newly created `data` folder, `tri-allelic` with `1000` haplotypes and `25000` SNPs; one with no missing data and a copy of the first one with 5\% missing data.

Example, insted:
```sh
./bin/gen data 2 1000 52000
```
Will generate a matrix under `data` folder, `bi-allelic` with `1000` haplotypes and `52000` SNPs.

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

