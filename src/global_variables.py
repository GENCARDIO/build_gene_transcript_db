import logging as log
import sys
import pyliftover
import os
import yaml

# Load the liftover chain file for GRCh38 to GRCh37
lo_38_to_37 = pyliftover.LiftOver('hg38', 'hg19')
lo_37_to_38 = pyliftover.LiftOver('hg19', 'hg38')
# Create a logging instance with a custom name
logging = log.getLogger("Mylogging")

# Create a formatter without the root prefix
log_formatter = log.Formatter('%(asctime)s - %(levelname)s - %(message)s')

# Create a StreamHandler and set its output to sys.stdout
console_handler = log.StreamHandler(sys.stdout)
console_handler.setLevel(log.DEBUG)
console_handler.setFormatter(log_formatter)

# Add the handler to the logging
logging.addHandler(console_handler)
logging.setLevel(log.DEBUG)

current_path = os.getcwd()

# Create a file handler:
file_handler = log.FileHandler("logfile.log")
file_handler.setLevel(log.DEBUG)
file_handler.setFormatter(log_formatter)

logging.addHandler(file_handler)

def get_conf(conf_file):
    with open(conf_file, "r") as file:
        data = yaml.safe_load(file)

    genes_db_name = data["database_name"]
    main_path = data["main_path"]
    refseq_gff_grch38 = os.path.join(
        main_path, data["resources"]["refseq"]["grch38"]["file_path"]
    )
    refseq_gff_grch37 = os.path.join(
        main_path, data["resources"]["refseq"]["grch37"]["file_path"]
    )
    mock_refseq_grch38 = os.path.join(
        main_path, data["resources"]["refseq"]["grch38"]["mock_file_path"]
    )
    mock_refseq_grch37 = os.path.join(
        main_path, data["resources"]["refseq"]["grch37"]["mock_file_path"]
    )
    mane_ensembl_gff = os.path.join(
        main_path, data["resources"]["mane"]["grch38"]["ensembl_file_path"]
    )
    mane_refseq_gff = os.path.join(
        main_path, data["resources"]["mane"]["grch38"]["refseq_file_path"]
    )
    mock_mane_ensembl_gff = os.path.join(
        main_path, data["resources"]["mane"]["grch38"]["mock_ensembl_file_path"]
    )
    mock_mane_refseq_gff = os.path.join(
        main_path, data["resources"]["mane"]["grch38"]["mock_refseq_file_path"]
    )
    gencode_gff_grch38 = os.path.join(
        main_path, data["resources"]["gencode"]["grch38"]["file_path"]
    )
    gencode_gff_grch37 = os.path.join(
        main_path, data["resources"]["gencode"]["grch37"]["file_path"]
    )
    mock_gencode_gff_grch38 = os.path.join(
        main_path, data["resources"]["gencode"]["grch38"]["mock_file_path"]
    )
    mock_gencode_gff_grch37 = os.path.join(
        main_path, data["resources"]["gencode"]["grch37"]["mock_file_path"]
    )
    lrg_gff = os.path.join(
        main_path, data["resources"]["lrg"]["file_path"]
    )
    hugo_path = os.path.join(
        main_path, data["resources"]["hugo"]["file_path"]
    )

    return (
        genes_db_name,
        refseq_gff_grch38,
        refseq_gff_grch37,
        mock_refseq_grch38,
        mock_refseq_grch37,
        mane_ensembl_gff,
        mane_refseq_gff,
        mock_mane_ensembl_gff,
        mock_mane_refseq_gff,
        gencode_gff_grch38,
        gencode_gff_grch37,
        mock_gencode_gff_grch38,
        mock_gencode_gff_grch37,
        lrg_gff,
        hugo_path,
        main_path
    )


def check_genome_version(version):
    """
    check genome version is 19, 37  or 38:

    Params:
    ------
        version: genome version
    Return:
    -------
        37(int): if version == 37 or 19
        38(int): if version == 38
        ValueError: If version is not 37, 38 or 19
    """
    allowed_versions = [19, 37, 38]
    if version not in allowed_versions:
        raise ValueError(
            f"Invalid genome version. Allowed versions are {allowed_versions}"
        )
    if version == 37 or version == 19:
        return 37
    elif version == 38:
        return (38)


def convert_chromosomes(chromosome: str):
    """
    convert Refseq chromosome nomenclature to a number.
    E.g. 'NC_000001' : '1'
    if chromosome not exists in chromosome_mapping keys, it
    returns the chromosome in refseq nomenclature
    
    Params:
    -------
        chromosome(str): chromosome name in refseq nomenclature
    
    Return:
    -------
        chromosome_mapping[chromosome] : chromosome number in string format
        chromosome: if chromosome not in chromosome mapping it returns the 
            refseq format
         """
    chromosome_mapping = {
        'NC_000001': 'chr1',
        'NC_000002': 'chr2',
        'NC_000003': 'chr3',
        'NC_000004': 'chr4',
        'NC_000005': 'chr5',
        'NC_000006': 'chr6',
        'NC_000007': 'chr7',
        'NC_000008': 'chr8',
        'NC_000009': 'chr9',
        'NC_000010': 'chr10',
        'NC_000011': 'chr11',
        'NC_000012': 'chr12',
        'NC_000013': 'chr13',
        'NC_000014': 'chr14',
        'NC_000015': 'chr15',
        'NC_000016': 'chr16',
        'NC_000017': 'chr17',
        'NC_000018': 'chr18',
        'NC_000019': 'chr19',
        'NC_000020': 'chr20',
        'NC_000021': 'chr21',
        'NC_000022': 'chr22',
        'NC_000023': 'chrX',
        'NC_000024': 'chrY',
        'NC_012920': 'chrM'
    }

    if chromosome in chromosome_mapping:
        return chromosome_mapping[chromosome]
    else:
        return chromosome