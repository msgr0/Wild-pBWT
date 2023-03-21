# `Wild-pBWT`

A PBWT-based algorithm for identifying all Maximal Perfect Haplotype Blocks with Wildcards (MPHBw).

## Build
SDSL by @SimonGog is required to build the project. More information available at [GitHub Repo](https://github.com/simongog/sdsl-lite)
Run
```
make
```
inside the main folder, to build `wild-pbtw`, `gen`, and `err`.

### `wild-pbwt`

The `wild-pbwt` binary under the `bin` subfolder computes all the maximal haplotype blocks with wildcards from a given input, with the extended pBWT

```sh
./bin/wild-pbwt -f <filename> -a <t-alleles> [-c|-o] [-v <verbose output>] [-b <block_size>] [-g <buffer_size>]
```
where:
- `<filename>` is the input file containing the haplotype panel in ASCII format where each line represents a variation site and each column is a single haplotype (i.e., similar to VCFs). Wildcards are represented with character `*`.
- `<t-alleles>` is the alphabet size (not counting `*`), i.e., the maximum number of alleles in a single site

If ran with `-o` flag, blocks will be output to standard output.
If ran with `-c` flag, blocks will not be output, only counted.
To specify a block_size, run with `-o` (`-c` not permitted)


### `gen`

The `gen` binary under the `bin` subfolder generates a matrix M x N with a specified wildcard data rate (if given)
```sh
./bin/gen <save_directory> <t-alleles> <haplotypes#> <SNPs#> <error_rate> 
```

## `err`

The `err` binary under the `bin` subfolder generates a matrix from an input, inserting a wildcard with probability <wild_rate>\% at each position.

```sh
./bin/gen <out_filename> <input_matrix> <wild_rate>
```


