import gzip
from src.global_variables import logging, convert_chromosomes


mane_gff = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/mane/MANE.GRCh38.v1.1.ensembl_genomic.gff.gz"
def get_mane(mane_gff: str):
    '''
        Function that reads a Mane gff3 file and returns two structured dictionaries
            with genes and transcripts
        :param str mane_gff: RefSeq gff3 file
        :return: two dicts, one with primary keys pointing to an enst_id, and a second one
            with gene id's as primary keys
        :rtype: Two defaultdict(dict)
    '''
    mane_select_ens_id = set()
    mane_clin_ens_id = set()
    mane_select_refseq_id = set()
    mane_clin_refseq_id = set()
    with gzip.open(mane_gff, 'rt') as fin:
        for line in fin:
            mane_clin = False
            mane_select = False

            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            feature = tmp[2].lower()
            refseq_chr = tmp[0]
            chr = convert_chromosomes(refseq_chr)
            pos = tmp[3]
            end = tmp[4]
            info= tmp[8].split(";")
            
            if feature == "transcript":
                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    if field == "id":
                        trans_id = result
                    if field == "tag":
                        tags = result.split(",")
                        for tag in tags:
                            if tag == "MANE_Select":
                                mane_select = True
                            if tag == "MANE_Plus_Clinical":
                                mane_clin = True
                    if field == "dbxref":
                        refseq_trans_id = result.split(":")[1]
                if mane_clin is True:
                    mane_clin_ens_id.add(trans_id)
                    mane_clin_refseq_id.add(refseq_trans_id)
                if mane_select is True:
                    mane_select_ens_id.add(trans_id)
                    mane_select_refseq_id.add(refseq_trans_id)
    
    return (
        mane_clin_ens_id,
        mane_select_ens_id,
        mane_clin_refseq_id,
        mane_select_refseq_id
    )


(
    mane_clin_ens_id,
    mane_select_ens_id,
    mane_clin_refseq_id,
    mane_select_refseq_id
) = get_mane(mane_gff)
