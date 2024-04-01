# m6Ad-SNVs

A Python pipeline for the analysis of variant-dependent m6A modifications within the Human genome.

## Requirements

This pipeline requires Python 3.8 with the following packages:

- _biopython_ (>=1.83)
- _pandas_ (>=2.2.1)
- _pysam_ (>=0.22.0)
- _tqdm_ (>=4.38.0)
- _wget_ (>=3.2)

They can be installed all at once by running `pip` over the requirements file under the `src` folder as follow: 

```bash
pip install src/requirements.txt
```

It also requires `RNAFold` as an external software dependency as part of the [ViennaRNA package](https://www.tbi.univie.ac.at/RNA/). It can be easily installed with `conda` by typing the following command in your terminal:

```bash
conda install viennarna
```

> __Note:__ 
> you need to add the `bioconda` channel to your conda instance before running installing the `viennarna` package. If you are wondering how to add the `bioconda` channel, have a look at the [Bioconda](https://bioconda.github.io) website.

## Usage

Once all the requirements are installed, you can finally run the pipeline which is all enclosed into the `m6Ad-SNVs.py` Python script. Here is a list of available options:

| Option             | Default            | Mandatory | Description  |
|:-------------------|:-------------------|:---------:|:-------------|
| `--bed`            |                    | ⚑         | Path to the BED file with the coordinates of the genomic regions of interest |
| `--bed-skip-lines` | `0`                |           | Skip the first number of lines in the input BED file |
| `--vcf`            |                    | ⚑         | Path to the VCF file with the SNPs coordinates |
| `--genome`         |                    | ⚑         | Path to the fasta genome file |
| `--search`         | `[AGT][AG]AC[ACT]` |           | Search for a specific pattern in the sequences. It accepts regular expressions. It searches for DRACH sites by default |
| `--out-table`      |                    |           | Path to the output file (optional, it overwrites the output file if it exists) |
| `--out-html`       |                    |           | Path to the output folder with RNAfold structures (optional, it overwrites the output files if they exist) |
| `--nproc`          | `1`                |           | Make it parallel |
| `--version`        |                    |           | Print the tool version and exit |

## Credits

Please credit our work in your manuscript by citing:

> _Manuscript in preparation_

## Support

If you need support, please open an [Issue](https://github.com/cumbof/m6Ad-SNVs/issues).