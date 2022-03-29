# GENE_TRANSCRIPT_DATABASE

## Installation
This code has been successfully tested to run in python 3.8.
You can check your python3 version by prompting: `python3 -V`.
Also, you will need `pip3` as a package manager.

It is advisable to use a virtual environment to build the project:
`python3 -m virtualenv mypython`

### 1. Get CrossMap

CrossMap is leveraged to perform coordinate conversions between genome versions.
To install it:

`pip3 install CrossMap`

You can check that CrossMap is installed properly by prompting on the terminal: `python3 CrossMap.py`

### 2. Install python3 dependencies
`pip3 install -r requeriments.txt`

## Basic usage:
Run the code with this command. You may want to introduce an output directory using parameter `--output_dir`.

`python3 gene_isoform_models.py --output_dir <output_directory>`

Upon execution, you will see this folder structure:

```
GENE_TRANSCRIPT_DATABASE
│
└───MANE
│   └───1.0
│       │───GRCh37
│             MANE.GRCh37.v1.0.ensembl_genomic.gff.gz
│       │───GRCh38
│             MANE.GRCh38.v1.0.ensembl_genomic.gff.gz
│
└───ENSEMBL
│   └───105
│       │───GRCh37
│             Homo_sapiens.GRCh37.105.gff3.gz
│       │───GRCh38
│             Homo_sapiens.GRCh38.105.gff3.gz
│
└───RefSeq
│   └───102
│       │───GRCh37
│             GRCh37_latest_genomic.gff.gz
│       │───GRCh38
│             GRCh38_latest_genomic.gff.gz
│
└───LRG
│   └───1.10
│       list_LRGs_transcripts_xrefs.txt
```
