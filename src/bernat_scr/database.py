import os
import sys
import gzip
import re
from collections import defaultdict
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker
from sqlalchemy_utils import database_exists, create_database
from sqlalchemy import Column, Integer, String, Float
from datetime import datetime
import logging
logger = logging.getLogger(__name__)

def get_mane_transcripts(mane_gff: str):
    '''
        Function that reads a MANE gff3 file and returns an structured dictionary
        :param str mane_gff: MANE gff3 file
        :return: A dict with primary keys pointing to an enst_id,
            and secondary keys with MANE information
        :rtype: A defaultdict(dict)
    '''
    mane_transcripts_dict = defaultdict(dict)
    with gzip.open(mane_gff, 'rt') as fin:
        for line in fin:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            feature = tmp[2]
            chr  = tmp[0]
            pos  = tmp[3]
            end  = tmp[4]
            info = tmp[8].split(";")
            if feature == "transcript":
                enst_id_full = next(filter(lambda item: item.startswith("transcript_id="), info), None)
                enst_id_full = enst_id_full.replace("transcript_id=", "")
                tmp_enst = enst_id_full.split(".")
                enst_id = tmp_enst[0]
                tag = next(filter(lambda item: item.startswith("tag="), info), None)
                tag = tag.replace("tag=", "")
                refseq_id = next(filter(lambda item: item.startswith("Dbxref=RefSeq:"), info), None)
                refseq_id = refseq_id.replace("Dbxref=RefSeq:", "")
                if not refseq_id:
                    refseq_id = "."
                if not enst_id in mane_transcripts_dict:
                    mane_transcripts_dict[enst_id] = defaultdict(dict)
                    mane_transcripts_dict[enst_id]['MANE_Select'] = "."
                    mane_transcripts_dict[enst_id]['MANE_Plus_Clinical'] = "."
                    if 'MANE_Select' in tag:
                        mane_transcripts_dict[enst_id]['MANE_Select'] = refseq_id
                    if 'MANE_Plus_Clinical' in tag:
                        mane_transcripts_dict[enst_id]['MANE_Plus_Clinical'] = refseq_id

    return mane_transcripts_dict

def get_gencode_transcripts(gencode_gff: str):
    '''
        Function that reads a GENCODE gff3 file and returns two structured dictionaries
            with genes and transcripts
        :param str mane_gff: GENCODE gff3 file
        :return: two dicts, one with primary keys pointing to an enst_id, and a second one
            with gene id's as primary keys
        :rtype: Two defaultdict(dict)
    '''
    gencode_genes = defaultdict(dict)
    gencode_genes_dict = defaultdict(dict)
    gencode_transcripts_dict = defaultdict(dict)

    with gzip.open(gencode_gff, 'rt') as fin:
        for line in fin:
            line = line.rstrip("\n")

            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            feature = tmp[2]
            chr  = tmp[0].replace("chr", "")
            pos  = tmp[3]
            end  = tmp[4]
            info = tmp[8].split(";")

            if feature == "CDS":
                gene_name = next(filter(lambda item: item.startswith("gene_name="), info), None)
                gene_name = gene_name.replace("gene_name=", "")
                enst_id = next(filter(lambda item: item.startswith("transcript_id="), info), None)
                if enst_id:
                    enst_id = enst_id.replace("transcript_id=", "")
                    tmp = enst_id.split(".")
                    enst_id = tmp[0]
                else:
                    enst_id = next(filter(lambda item: item.startswith("Parent=transcript:"), info), None)
                    enst_id = enst_id.replace("Parent=transcript:", "")
                    tmp = enst_id.split(".")
                    enst_id = tmp[0]
                ensp_id = next(filter(lambda item: item.startswith("protein_id="), info), None)
                ensp_id = ensp_id.replace("protein_id=", "")

                if not gene_name in gencode_genes_dict:
                    gencode_genes_dict[gene_name]['transcripts'] = set()
                gencode_genes_dict[gene_name]['transcripts'].add(enst_id)

                if not enst_id in gencode_transcripts_dict:
                    gencode_transcripts_dict[enst_id] = defaultdict(dict)
                    gencode_transcripts_dict[enst_id]['exons'] = defaultdict(dict)
                coordinates = ("{}:{}-{}").format(chr, pos, end)
                gencode_transcripts_dict[enst_id]['exons'][coordinates] = ""
                gencode_transcripts_dict[enst_id]['gene_name'] = gene_name
                gencode_transcripts_dict[enst_id]['ensp'] = ensp_id

    return gencode_transcripts_dict, gencode_genes_dict

def get_ensembl_transcripts(ensembl_gff: str):
    '''
    '''
    ensembl_genes = defaultdict(dict)
    ensembl_genes_dict = defaultdict(dict)

    with gzip.open(ensembl_gff, 'rt') as fin:
        for line in fin:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            feature = tmp[2]
            chr  = tmp[0]
            pos  = tmp[3]
            end  = tmp[4]
            info = tmp[8].split(";")
            if feature == "mRNA":
                gene_name = next(filter(lambda item: item.startswith("Name="), info), None)
                if gene_name:
                    tmp_gene_name = gene_name.split("-")
                    gene_name = tmp_gene_name[0].replace("Name=", "")
                else:
                    gene_id = next(filter(lambda item: item.startswith("Parent=gene:"), info), None)
                    gene_id = gene_id.replace("Parent=gene:", "")
                    gene_name = gene_id
                enst_id = next(filter(lambda item: item.startswith("transcript_id="), info), None)
                enst_id = enst_id.replace("transcript_id=", "")

                if not enst_id in ensembl_genes:
                    ensembl_genes[enst_id] = gene_name
                if not gene_name in ensembl_genes_dict:
                    ensembl_genes_dict[gene_name] = defaultdict(dict)
                    ensembl_genes_dict[gene_name]['ENSP'] = ""
                    ensembl_genes_dict[gene_name]['transcripts'] = []
                ensembl_genes_dict[gene_name]['transcripts'].append(enst_id)

    ensembl_transcripts = defaultdict(dict)
    with gzip.open(ensembl_gff, 'rt') as fin:
        for line in fin:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            feature = tmp[2]
            chr  = tmp[0]
            pos  = tmp[3]
            end  = tmp[4]
            info = tmp[8].split(";")

            # Protein coding transcript
            if feature == "CDS":
                enst = next(filter(lambda item: item.startswith("Parent=transcript:"), info), None)
                enst = enst.replace("Parent=transcript:", "")

                ensp = next(filter(lambda item: item.startswith("protein_id="), info), None)
                ensp = ensp.replace("protein_id=", "")

                if not enst in ensembl_transcripts:
                    ensembl_transcripts[enst] = defaultdict(dict)
                    ensembl_transcripts[enst]['exons'] = defaultdict(dict)
                coordinates = ("{}:{}-{}").format(chr, pos, end)
                ensembl_transcripts[enst]['exons'][coordinates] = ""
                ensembl_transcripts[enst]['gene_name'] = ensembl_genes[enst]
                ensembl_transcripts[enst]['ensp'] = ensp

    return ensembl_transcripts, ensembl_genes_dict

def get_refseq_transcripts(refseq_gff: str):
    '''
        Function that reads a RefSeq gff3 file and returns two structured dictionaries
            with genes and transcripts
        :param str mane_gff: RefSeq gff3 file
        :return: two dicts, one with primary keys pointing to an enst_id, and a second one
            with gene id's as primary keys
        :rtype: Two defaultdict(dict)
    '''
    refseq_transcripts = defaultdict(dict)
    refseq_genes_dict  = defaultdict(dict)
    with gzip.open(refseq_gff, 'rt') as fin:
        for line in fin:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            feature = tmp[2]
            chr = tmp[0]
            pos = tmp[3]
            end = tmp[4]
            info= tmp[8].split(";")
            if feature == "CDS":
                gene_name = next(filter(lambda item: item.startswith("gene="), info), None)
                isoform  = next(filter(lambda item: item.startswith("Parent="), info), None)
                isoform = isoform.replace("Parent=rna-", "")

                protein_id = next(filter(lambda item: item.startswith("protein_id="), info), None)
                if protein_id:
                    protein_id = protein_id.replace("protein_id=", "")
                else:
                    protein_id = "."

                if gene_name:
                    gene_name = gene_name.replace("gene=", "")
                else:
                    gene_name = isoform
                if not isoform in refseq_transcripts:
                    refseq_transcripts[isoform] = defaultdict(dict)
                    refseq_transcripts[isoform]['exons'] = defaultdict(dict)

                coordinates = ("{}:{}-{}").format(chr, pos, end)
                refseq_transcripts[isoform]['exons'][coordinates] = ""
                refseq_transcripts[isoform]['gene_name'] = gene_name
                refseq_transcripts[isoform]['protein_id']= protein_id
                if not gene_name in refseq_genes_dict:
                    refseq_genes_dict[gene_name] = defaultdict(dict)
                    refseq_genes_dict[gene_name]['transcripts'] = set()
                refseq_genes_dict[gene_name]['transcripts'].add(isoform)

    return refseq_transcripts, refseq_genes_dict

def get_lrg_transcripts(lrg_txt: str):
    '''
        Function that reads a plain txt file with LRG gene/transcript id's gff3 and
            returns a dictioanry
        :param str lrg_txt: LRG txt file
        :return: A dictcionary with primary keys pointing to an enst_id and secondary
            keys with LRG id's
        :rtype: defaultdict(dict)
    '''
    lrg_transcripts = defaultdict(dict)
    with open(lrg_txt) as f:
        for line in f:
        # Last modified: 30-03-2021@22:00:06
        # LRG	HGNC_SYMBOL	REFSEQ_GENOMIC	LRG_TRANSCRIPT	REFSEQ_TRANSCRIPT	ENSEMBL_TRANSCRIPT	CCDS
        # LRG_1	COL1A1	NG_007400.1	t1	NM_000088.3	ENST00000225964.10	CCDS11561.1
            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            lrg_id = tmp[0]
            lrg_transcript = tmp[3]
            tmp_refseq = tmp[4].split(".")
            refseq_id  = tmp_refseq[0]
            tmp_enst   = tmp[5].split(".")
            enst_id    = tmp_enst[0]
            if not enst_id in lrg_transcripts:
                lrg_transcripts[enst_id] = defaultdict(dict)
                lrg_transcripts[enst_id]['LRG_ID'] = lrg_id
                lrg_transcripts[enst_id]['LRG_TRANSCRIPT'] = lrg_transcript

        f.close()
    return lrg_transcripts

def build_database(genome_version: str, input_file_dict, output_dir, versions_dict):
    '''
        Build a sqlite database with gene/transcript identifiers from GENCODE, RefSeq,
        MANE and LRG.
        :param str genome_version: human genome version
        :param dict input_file_dict: dict with gff3 file location
        :param str output_dir: output directory
        :param dict versions_dict: dict with resource versions
    '''

    mane_transcripts_dict = get_mane_transcripts(input_file_dict[genome_version]['mane'])
    gencode_transcripts_dict, gencode_genes_dict \
        = get_gencode_transcripts(input_file_dict[genome_version]['gencode'])
    refseq_transcripts_dict, refseq_genes_dict   \
        = get_refseq_transcripts(input_file_dict[genome_version]['refseq'])
    lrg_transcripts = get_lrg_transcripts(input_file_dict[genome_version]['lrg'])

    missing_genes =  {
        'refseq'  : [],
        'ensembl' : [],
        'gencode' : []
    }

    # Now define and create the genes database
    genes_db= output_dir + "/genes.db"
    engine  = create_engine('sqlite:///' + genes_db)
    Session = sessionmaker(bind=engine)
    session = Session()
    Base = declarative_base()

    class TranscriptModels(Base):
        __tablename__ = 'TRANSCRIPT_MODELS'
        ID = Column(Integer, primary_key=True)
        GENE_SYMBOL = Column(String, nullable=False)
        ENST_ID = Column(String)
        ENSP_ID = Column(String)
        REFSEQ_CDS_COMPLETE = Column(String)
        REFSEQ_PROTEIN_ID = Column(String)
        MANE_SELECT = Column(String)
        MANE_PLUS_CLINICAL = Column(String)
        LRG_ID = Column(String)
        LRG_TRANSCRIPT= Column(String)

    class Releases(Base):
        __tablename__ = 'RELEASES'
        ID = Column(Integer, primary_key=True)
        DATE = Column(String)
        GENCODE = Column(String)
        REFSEQ = Column(String)
        MANE = Column(String)
        LRG = Column(String)

    Base.metadata.create_all(engine)

    now = datetime.now() # current date and time
    date_time = now.strftime("%d/%m/%Y")

    # add source versions
    releases = Releases(DATE=date_time,
        GENCODE=versions_dict[genome_version]['gencode'],
        # ENSEMBL=versions_dict[genome_version]['ensembl'],
        REFSEQ=versions_dict[genome_version]['refseq'],
        MANE=versions_dict[genome_version]['mane'],
        LRG=versions_dict[genome_version]['lrg'])
    session.add(releases)
    session.commit()

    written_records =  0
    objects = []
    for gene in gencode_genes_dict:
        for enst_id in gencode_genes_dict[gene]['transcripts']:
            tag = "."
            mane_select = "."
            mane_plus_clinical = "."
            ensp = "."
            refseq_protein_id = "."
            if enst_id in gencode_transcripts_dict:
                if 'ensp' in gencode_transcripts_dict[enst_id]:
                    ensp = gencode_transcripts_dict[enst_id]['ensp']
            refseq_protein_id = "."
            lrg_id = "."
            lrg_transcript = "."
            if enst_id in lrg_transcripts:
                lrg_id = lrg_transcripts[enst_id]['LRG_ID']
                lrg_transcript = lrg_transcripts[enst_id]['LRG_TRANSCRIPT']
            if enst_id in mane_transcripts_dict:
                mane_select = mane_transcripts_dict[enst_id]['MANE_Select']
                mane_plus_clinical = mane_transcripts_dict[enst_id]['MANE_Plus_Clinical']

            if not 'exons' in gencode_transcripts_dict[enst_id]:
                missing_genes['gencode'].append(gene)
                continue

            enst_exons = gencode_transcripts_dict[enst_id]['exons']
            if not gene in refseq_genes_dict:
                missing_genes['refseq'].append(gene)
                continue

            identical_refseq = []
            for refseq_id in refseq_genes_dict[gene]['transcripts']:
                refseq_exons = refseq_transcripts_dict[refseq_id]['exons']
                refseq_protein_id = refseq_transcripts_dict[refseq_id]['protein_id']
                if enst_exons == refseq_exons:
                    identical_refseq.append(refseq_id)

            transcript_model = TranscriptModels(GENE_SYMBOL=gene,
                ENST_ID=enst_id,
                ENSP_ID=ensp,
                REFSEQ_CDS_COMPLETE=','.join(identical_refseq),
                REFSEQ_PROTEIN_ID=refseq_protein_id,
                MANE_SELECT=mane_select,
                MANE_PLUS_CLINICAL=mane_plus_clinical,
                LRG_ID=lrg_id,
                LRG_TRANSCRIPT=lrg_transcript)
            try:
                objects.append(transcript_model)
            except:
                msg = "ERROR: An error occured at database writing"
                logging.error(msg)
                raise SQLAlchemyError(msg)
            else:
                written_records+=1
    try:
        session.bulk_save_objects(objects)
        session.commit()
    except:
        msg = "ERROR: An error occured at database writing"
        logging.error(msg)
        raise SQLAlchemyError(msg)
    else:
        msg = "INFO: Written up to {} records successfully".format(written_records)
        logging.info(msg)
