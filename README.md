# build_gene_transcript_db

Create a database to link gene names (HUGO) with the following genome resources: GENCODE/Ensembl, RefSeq, MANE and LRG.

## TODO
Check disparity between transcript versions (i.e ENST00001.X)

## Installation
This code has been successfully tested to run in python 3.8.
You can check your python3 version by prompting: `python3 -V`.
Also, you will need `pip3` as a package manager.

It is advisable to use a virtual environment to build the project:
`python3 -m virtualenv mypython`

### 1. Get CrossMap

CrossMap handles coordinate conversions between genome versions.
To install it:

`pip3 install CrossMap`

You can check that CrossMap is installed properly by prompting on the terminal: `CrossMap.py`

### 2. Install python3 dependencies
`pip3 install -r requeriments.txt`

## Basic usage:
Run the code with this command. You may want to introduce an output directory using parameter `--output_dir`.

`python3 gene_isoform_models.py --output_dir <output_directory>`

Upon execution, the following database will be generated:

| **GENE_SYMBOL** | **ENST_ID**     | **ENSP_ID**       | **REFSEQ_CDS_COMPLETE** | **REFSEQ_PROTEIN_ID** | **MANE_SELECT** | **MANE_PLUS_CLINICAL** | **LRG_ID** | **LRG_TRANSCRIPT** |
|-----------------|-----------------|-------------------|-------------------------|-----------------------|-----------------|------------------------|------------|--------------------|
| **SCN5A**       | ENST00000423572 | ENSP00000398266.2 | NM_000335.5             | NP_000326.2           | NM_000335.5     | .                      | LRG_289    | t2                 |
| **SCN5A**       | ENST00000327956 | ENSP00000333674.6 |                         | NP_000326.2           | .               | .                      | .          | .                  |
| **SCN5A**       | ENST00000414099 | ENSP00000398962.2 | NM_001099405.2          | NP_000326.2           | .               | .                      | .          | .                  |
| **SCN5A**       | ENST00000333535 | ENSP00000328968.4 | NM_198056.3             | NP_000326.2           | .               | .                      | LRG_289    | t1                 |
| **SCN5A**       | ENST00000413689 | ENSP00000410257.1 | NM_001099404.2          | NP_000326.2           | .               | NM_001099404.2         | LRG_289    | t3                 |
| **SCN5A**       | ENST00000450102 | ENSP00000403355.2 | NM_001160161.2          | NP_000326.2           | .               | .                      | .          | .                  |
| **SCN5A**       | ENST00000449557 | ENSP00000413996.2 |                         | NP_000326.2           | .               | .                      | .          | .                  |
| **SCN5A**       | ENST00000455624 | ENSP00000399524.2 | NM_001160160.2          | NP_000326.2           | .               | .                      | .          | .                  |


All gff3 files that populate the database will be downloaded within the following folder schema:

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
└───GENCODE
│   └───39
│       │───GRCh37
│             gencode.v39lift37.annotation.gff3.gz
│       │───GRCh38
│             gencode.v39.annotation.gff3.gz
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
