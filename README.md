# <span style="color:red;">sp</span><span style="color:orange;">ee</span><span style="color:yellow;">dd</span><span style="color:lime;">em</span><span style="color:blue;">O</span><span style="color:indigo;">N</span><span style="color:violet;">T</span>

The **spee**dy **dd**-RADseq **dem**ultiplexer for **ONT** (and Illumina) data.

## Purpose

A `Python` program designed to demultiplex ONT-barcoded ddRADseq data.
This is the sister program which also demuxxes Illumina data.
Currently being developed for [Oikos](https://oikosgenomics.org)'s genomics low-cost parentage analyses.

## Quickstart

The minimal example below will write all demuxxed fastq files to `speeddemONT_out` in your working directory.

```
speeddemONT -f <input.fq.fz> -d <input_demux_construct_file.tsv> -fa .9 -b 9
```

### Options

| Option | Default | Data type | Description |
| -- | -- | -- | -- |
| `-h`, `--help` | `FALSE` | Flag | Print a help message and exit. |
| `-f`, `--fastq` | `null` | String | Path to the input FASTQ file to be demultiplexed. Mandatory. |
| `-d`, `--demux` | `null` | String | Path to the demux construct TSV file. Mandatory. |
| `-b`, `--buffer` | `0` | Int | The integer number of base pairs outside of the long element the short element can align to. Defaults to 0 (internal matching only). Mandatory. |
| `-p`, `--prefix` | `speeddemONT_out` | String | Output directory prefix. Demultiplexed FASTQ files will be written to this directory. |
| `-fa`, `--fuzzy_aln_percent` | `0.9` | Float | The minimum percent identity (0.0-1.0) needed to fuzzy-match a full index or barcode to a sequence. Used for `index_full` and `barcode_full` alignments. |
| `-ea`, `--exact_aln_percent` | `1.0` | Float | The minimum percent identity (0.0-1.0) needed to exact-match a short index or barcode to a sequence. Used for `index` and `barcode` alignments. |
| `-ea`, `--exact_aln_percent` | `1.0` | Float | The minimum percent identity (0.0-1.0) needed to exact-match a short index or barcode to a sequence. Used for `index` and `barcode` alignments. |


## Inputs
