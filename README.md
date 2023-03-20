# `Wild-pBWT`

A PBWT-based algorithm for identifying all Maximal Perfect Haplotype Blocks with Wildcards (MPHBw).

## Build
SDSL by @SimonGog is required to build the project. More information available at [GitHub Repo](https://github.com/simongog/sdsl-lite)
Run
```
make
```
inside the main folder, to build `pbtw`, `gen`, and `err`.

### `pbwt`
Run
```
make pbwt
```
to build the executable for `pbwt` only.

The `pbwt` binary under the `bin` subfolder computes all the MPHBw from a given input, with the extended pBWT
```sh
./bin/gen -f <filename> -a <t-alleles> -c <y:count blocks> -v <verbose output> -o <y:compute blocks> -b <block_size> -g <buffer_size>
```
If ran with `-o` flag, blocks will be outputted to std:cout.
If ran with `-c` flag, blocks wont be outputted, only counted.
To specify a block_size, run with -o (-c not permitted)


### `gen`
Run
```
make gen
```
to build the executable for `gen` only.
The `gen` binary under the `bin` subfolder generates a matrix M x N with a specified missing data rate (if given)
```sh
./bin/gen <save_directory> <t-alleles> <haplotypes#> <SNPs#> <error_rate> 
```

## 'err`
Run
```
make err
```
to build the executable for `err` only.
The `err` binary under the `bin` subfolder generates a matrix from an input, inserting <error_rate>\% wildcards
```sh
./bin/gen <out_filename> <input_matrix> <error_rate> 
```


