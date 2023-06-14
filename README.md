# m6A

The `m6A.py` script requires one input one. It currently uses the [m6A table](https://rmvar.renlab.org/download/RMVar_Human_basic_info_m6A.txt) from the RMVar database. For each of the m6A entry (row) in the input table, RMVar reports a reference and alterative sequence.

> **Note:**
> It can automatically retrieve the input table by specifying the `--download` argument.
>
> More in general, it should be able to consider a bed file with a set of genomic sites and a VCF file with a series of SNPs.

The script integrates [RNAfold](), which is part of the [ViennaRNA](https://github.com/ViennaRNA/ViennaRNA) suite, for predicting the secondary structures of both the reference and the alterative sequences.

It currently considers m6A entries for which no SNP modifications occur on any DRACH sites on the reference (and alterative) sequence. The m6A must also occur on exons.

> **Note:**
> DRACH sites are detected by searching on the reference sequence with the following regular expression: `[AGT][AG]AC[ACT]`.
>
> Given a specific sequence, the script also implement the `--expand-left` and `--expand-righ` arguments for automatically retrieving the gene sequence through the NCBI Datasets utility and expand the reference (and implicitely the alterative) sequence _N_ bases to the left and _N_ bases to the right.

It finally compares the reference and alterative structures in dot-bracket notation to check whether the SNP modification changed the status of at least one DRACH site from paired to unpaired or the other way around.

Results are reported in a tab-delimited table by specifying `--table`. The script can also produce an HTML page for each of the target m6A sites with the graphical representation of both the reference and alterative structures with the [fornac](https://github.com/ViennaRNA/fornac) Javascript library by specifying `--structures`.

Here is the list of arguments that can also be inspected by typing `python m6A.py --help`:

| Argument            | Description  |
|:--------------------|:-------------|
| `--table`           | Path to the input table with m6A sites from the RMVar database |
| `--download`        | Automatically downlowd the m6A table from RMVar if `--table` is not specified |
| `--expand-left`     | Automatically retrieve the gene sequence from NCBI Datasets and expand the reference sequence by _N_ bases to the left |
| `--expand-right`    | Same as `--expand-left` but expands the sequence by _N_ bases to the right |
| `--paired-unpaired` | Report the m6A sites for which the SNP modifications caused a disruption of one of the DRACH sites in the reference structure (from paired to unpaired) |
| `--unpaired-paired` | Same as `--paired-unpaired` but considering an unpaired-to-paired status change of one of the DRACH sites in the reference sequence |
| `--out-table`       | Path to the output tab-delimited table with target m6A sites |
| `--out-structures`  | Path to the folder with an HTML page for each target m6A site with the graphical representation of both the reference and alterative structures |
| `--version`         | Print the script version on the stdout and exit |