import gzip
import os
import re
from src.global_variables import logging, convert_chromosomes

# file paths
mane_ensembl_gff = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/mane/MANE.GRCh38.v1.2.ensembl_genomic.gff.gz"
mane_refseq_gff = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/mane/MANE.GRCh38.v1.2.refseq_genomic.gff.gz"
mock_mane_ensembl_gff = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/mane/mock.MANE.GRCh38.v1.2.ensembl_genomic.gff.gz"
mock_mane_refseq_gff = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/mane/mock.MANE.GRCh38.v1.2.refseq_genomic.gff.gz"


def get_mane_conf(mane_gff_path, versions_dict, genome_version):
    """
    Extract mane versions of the files used to create the database

    Params:
    -------
        mane_gff_path: path of the mane file
        versions_dict: dict to store resources information
        genome_version: genome version (38 or 37)
    Return:
    -------
        versions_dict: resource information with mane info
    """
    version_pattern = r'\.v(\d+\.\d+)\.'
    filename = os.path.basename(mane_gff_path)

    match = re.search(version_pattern, filename)
    if "resources" not in versions_dict:
        versions_dict["resources"] = {}
    if "mane" not in versions_dict["resources"]:
        versions_dict["resources"]["mane"] = {}
    if "grch38" not in versions_dict["resources"]["mane"]:
        versions_dict["resources"]["mane"]["grch38"] = {}
    if "refseq" not in versions_dict["resources"]["mane"]["grch38"]:
        versions_dict["resources"]["mane"]["grch38"]["refseq"] = {}
    if "ensembl" not in versions_dict["resources"]["mane"]["grch38"]:
        versions_dict["resources"]["mane"]["grch38"]["ensembl"] = {}
    if "grch37" not in versions_dict["resources"]["mane"]:
        versions_dict["resources"]["mane"]["grch37"] = {}
    if "refseq" not in versions_dict["resources"]["mane"]["grch37"]:
        versions_dict["resources"]["mane"]["grch37"]["refseq"] = {}
    if "ensembl" not in versions_dict["resources"]["mane"]["grch37"]:
        versions_dict["resources"]["mane"]["grch37"]["ensembl"] = {}

    if match:
        version = match.group(1)
    else:
        version = "not detected"
    if "refseq" in filename:
        if genome_version == 37:
            versions_dict["resources"]["mane"]["grch37"]["refseq"]["version"] = version
            versions_dict["resources"]["mane"]["grch37"]["refseq"]["file_path"] = mane_gff_path
        elif genome_version == 38:
            versions_dict["resources"]["mane"]["grch38"]["refseq"]["version"] = version
            versions_dict["resources"]["mane"]["grch38"]["refseq"]["file_path"] = mane_gff_path
    elif "ensembl" in filename:
        if genome_version == 37:
            versions_dict["resources"]["mane"]["grch37"]["ensembl"]["version"] = version
            versions_dict["resources"]["mane"]["grch37"]["ensembl"]["file_path"] = mane_gff_path
        elif genome_version == 38:
            versions_dict["resources"]["mane"]["grch38"]["ensembl"]["version"] = version
            versions_dict["resources"]["mane"]["grch38"]["ensembl"]["file_path"] = mane_gff_path 

    return (versions_dict)

def get_ensembl_mane(mane_gff: str):
    '''
        Function that reads a Mane gff3 file and returns two structured dictionaries
            with genes and transcripts
        :param str mane_gff: RefSeq gff3 file
        :return: two dicts, one with primary keys pointing to an enst_id, and a second one
            with gene id's as primary keys
        :rtype: Two defaultdict(dict)
    '''

    mane_clin_ensembl_id = set()
    mane_select_refseq_ensembl = dict()
    mane_clin_refseq_id = set()
    mane_select_ensembl_refseq = dict()
    transcript_id_set = set()
    with gzip.open(mane_gff, 'rt') as fin:
        for line in fin:
            mane_clin = False
            mane_select = False

            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            feature = tmp[2].lower()
            info = tmp[8].split(";")
            
            if feature == "transcript":
                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1].upper()
                    if field == "id":
                        ensembl_trans_id = result
                        transcript_id_set.add(ensembl_trans_id)
                    if field == "tag":
                        tags = result.split(",")
                        for tag in tags:
                            if tag == "MANE_SELECT":
                                mane_select = True
                            if tag == "MANE_PLUS_CLINICAL":
                                mane_clin = True
                    if field == "dbxref":
                        if "-" in result:
                            result = result.split("-")[0]
                        refseq_trans_id = result.split(":")[1]

                if mane_clin is True:
                    mane_clin_ensembl_id.add(ensembl_trans_id)
                    mane_clin_refseq_id.add(refseq_trans_id)
                if mane_select is True:
                    mane_select_ensembl_refseq[ensembl_trans_id] = refseq_trans_id
                    mane_select_refseq_ensembl[refseq_trans_id] = ensembl_trans_id
    return (
        mane_clin_ensembl_id,
        mane_clin_refseq_id,
        mane_select_ensembl_refseq,
        mane_select_refseq_ensembl
    )


def get_refseq_mane(
        mane_gff: str,
        mane_clin_refseq_id: set,
        mane_clin_ensembl_id: set,
        mane_select_refseq_ensembl: dict,
        mane_select_ensembl_refseq: dict
):
    '''
        Function that reads a Mane gff3 file and returns two structured dictionaries
            with genes and transcripts
        :param str mane_gff: RefSeq gff3 file
        :return: two dicts, one with primary keys pointing to an enst_id, and a second one
            with gene id's as primary keys
        :rtype: Two defaultdict(dict)
    '''
    feature_set = set()
    with gzip.open(mane_gff, 'rt') as fin:
        for line in fin:
            mane_clin = False
            mane_select = False

            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            feature = tmp[2].lower()
            info = tmp[8].split(";")
            feature_set.add(feature)
            if feature == "mrna":
                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1].upper()
                    if field == "id":
                        if "RNA-" in result:
                            result = result.replace("RNA-", "")
                            if "-" in result:
                                result = result.split("-")[0]
                        refseq_trans_id = result
                    if field == "tag":
                        tags = result.split(",")  # TAGS ARE IN UPPERCASE
                        for tag in tags:
                            tag = tag.replace(" ", "_")
                            if "MANE_SELECT" in tag:
                                mane_select = True
                            if "MANE_PLUS_CLINICAL" in tag:
                                mane_clin = True
                    if field == "dbxref":
                        result = result.split(",")
                        for item in result:
                            if "ENSEMBL" in item:
                                ensembl_trans_id = item.split(":")[1]
                if mane_clin is True:
                    mane_clin_ensembl_id.add(ensembl_trans_id)
                    mane_clin_refseq_id.add(refseq_trans_id)
                if mane_select is True:
                    if ensembl_trans_id in mane_select_ensembl_refseq:
                        refseq_assigned = mane_select_ensembl_refseq[ensembl_trans_id]
                        if refseq_assigned != refseq_trans_id:
                            logging.critical(
                                f"MANE refseq id assigned for ensembl_id {ensembl_trans_id}\
                                not match between mane refseq and mane ensembl: \
                                    {refseq_assigned} {refseq_trans_id}"
                            )
                    else:
                        mane_select_ensembl_refseq[ensembl_trans_id] = refseq_trans_id
                    if refseq_trans_id in mane_select_refseq_ensembl:
                        ensembl_assigned = mane_select_refseq_ensembl[refseq_trans_id]
                        if ensembl_assigned != ensembl_trans_id:
                            logging.critical(
                                f"MANE ensembl id assigned for refseq_id {refseq_trans_id}\
                                not match between mane refseq and mane ensembl files: \
                                    {ensembl_assigned} {ensembl_trans_id}"
                            )
                    else:
                        mane_select_refseq_ensembl[refseq_trans_id] = ensembl_trans_id
        if len(mane_select_refseq_ensembl) != len(mane_select_ensembl_refseq):
            logging.critical(
                f"length of ensembl manes and refseq mane does not match:\
                ensembl: {len(mane_select_ensembl_refseq)}\
                refseq: {len(mane_select_refseq_ensembl)}"
            )
    logging.info(
        f"{len(mane_select_refseq_ensembl)} refseq ids and ensembl ids \
        extracted from {mane_ensembl_gff} and {mane_refseq_gff}")
    return (
        mane_clin_ensembl_id,
        mane_clin_refseq_id,
        mane_select_ensembl_refseq,
        mane_select_refseq_ensembl
    )


if "__main__" == __name__:

    (
        mane_clin_ensembl_id,
        mane_select_ens_id,
        mane_clin_refseq_id,
        mane_select_refseq_id
    ) = get_ensembl_mane(mane_ensembl_gff)

    (
        mane_clin_ensembl_id,
        mane_select_ens_id,
        mane_clin_refseq_id,
        mane_select_refseq_id
    ) = get_refseq_mane(
        mane_refseq_gff,
        mane_clin_refseq_id,
        mane_select_refseq_id
    )