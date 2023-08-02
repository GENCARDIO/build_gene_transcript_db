import logging as log
import sys
import pyliftover

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

def check_genome_version(version):
    allowed_versions = [19, 37, 38]
    if version not in allowed_versions:
        raise ValueError(
            f"Invalid genome version. Allowed versions are {allowed_versions}"
        )
    if version == 37 or version == 19:
        return 37
    elif version == 38:
        return(38)
def convert_chromosomes(chromosome):
    chromosome_mapping = {
        'NC_000001': '1',
        'NC_000002': '2',
        'NC_000003': '3',
        'NC_000004': '4',
        'NC_000005': '5',
        'NC_000006': '6',
        'NC_000007': '7',
        'NC_000008': '8',
        'NC_000009': '9',
        'NC_000010': '10',
        'NC_000011': '11',
        'NC_000012': '12',
        'NC_000013': '13',
        'NC_000014': '14',
        'NC_000015': '15',
        'NC_000016': '16',
        'NC_000017': '17',
        'NC_000018': '18',
        'NC_000019': '19',
        'NC_000020': '20',
        'NC_000021': '21',
        'NC_000022': '22',
        'NC_000023': 'X',
        'NC_000024': 'Y',
        'NC_012920': 'M'
    }

    if chromosome in chromosome_mapping:
        return chromosome_mapping[chromosome]
    else:
        return chromosome