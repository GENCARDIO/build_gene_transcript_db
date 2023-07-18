from dataclasses import dataclass, field
import os
import sys
import gzip
from collections import defaultdict
from collections import defaultdict
from sqlalchemy import create_engine, exc, exists
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.inspection import inspect
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy_utils import database_exists, create_database
from sqlalchemy import Column, Integer, String, Float, ForeignKey, Boolean, Table
import logging
import json

@dataclass
class EnsemblGenes:
    gene_id                 : str
    start                   : int
    end                     : int
    chromosome              : str
    id                      : str
    gene_name               : str
    gene_type               : str
    level                   : int
    tag                     : str = field(default=None)
    havana_gene             : str = field(default=None)
    hgnc_id                 : int = field(default=None)
    num_gencode_transcripts : int = field(default=0)
    num_refseq_transcripts  : int = field(default=0)
    refseq_same_gene_coords : str = field(default=None)

    def increment_num_gencode_transcripts(self):
        self.num_gencode_transcripts += 1
    
    def increment_num_refseq_transcripts(self):
        self.num_refseq_transcripts += 1
    
    def get_coords(self):
        return f"{self.chromosome}:{self.start}-{self.end}"


@dataclass
class RefseqGene:
    gene_id         : str
    start           : int
    end             : int
    chromosome      : int
    id              : str
    hgnc_id         : int
    name            : str
    gene_biotype    : str

    def get_coords(self):
        return f"{self.chromosome}:{self.start}-{self.end}"

@dataclass
class EnsemblTranscript:
    transcript_id       : str
    start               : int
    end                 : int
    chromosome          : str
    transcript_type     : str
    transcript_name     : str
    id                  : str
    gene_id             : str
    protein_id          : str = field(default=None)
    lrg_transcript      : str = field(default=None)
    lrg_id              : str = field(default=None)
    ccds                : str = field(default=None)
    numb_exons          : int = field(default=0)
    numb_cds            : int = field(default=0)
    mane_select         : bool = field(default=False)
    mane_clin           : bool = field(default=False)
    havana_transcript   : str = field(default=None)

    def increment_numb_exons(self):
        self.numb_exons += 1
    
    def increment_numb_cds(self):
        self.numb_cds += 1
    
    def get_coords(self):
        return f"{self.chromosome}:{self.start}-{self.end}"

@dataclass
class EnsemblExons:
    exon_id         : str
    id              : str
    exon_number     : int
    start           : int
    end             : int
    chromosome      : str
    transcript_id   : str

    def get_coords(self):
        return f"{self.chromosome}:{self.start}-{self.end}"

@dataclass
class EnsemblCds:
    id              : str
    start           : int
    end             : int
    chromosome      : str
    exon_id         : str
    transcript_id   : str
    exon_number     : int

    def get_coords(self):
        return f"{self.chromosome}:{self.start}-{self.end}"



# Configure the root logger
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Create a StreamHandler and set its output to sys.stdout
console_handler = logging.StreamHandler(sys.stdout)

# Optionally, you can set the logging level for the handler
console_handler.setLevel(logging.DEBUG)

# Add the handler to the root logger
logging.getLogger().addHandler(console_handler)

gene_table_cols = [
    "name",
    "id",
    "chr",
    "start",
    "end",
    "hgnc",
    "type"
]

gencode_gff = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/gencode/gencode.v43.annotation.gff3.gz"

def convert_chromosomes(chromosome):
    chromosome_mapping = {
        'NC_000001.11': '1',
        'NC_000002.12': '2',
        'NC_000003.12': '3',
        'NC_000004.12': '4',
        'NC_000005.10': '5',
        'NC_000006.12': '6',
        'NC_000007.14': '7',
        'NC_000008.11': '8',
        'NC_000009.12': '9',
        'NC_000010.11': '10',
        'NC_000011.10': '11',
        'NC_000012.12': '12',
        'NC_000013.11': '13',
        'NC_000014.9': '14',
        'NC_000015.10': '15',
        'NC_000016.10': '16',
        'NC_000017.11': '17',
        'NC_000018.10': '18',
        'NC_000019.10': '19',
        'NC_000020.11': '20',
        'NC_000021.9': '21',
        'NC_000022.11': '22',
        'NC_000023.11': 'X',
        'NC_000024.10': 'Y',
        'NC_012920.1': 'M'
    }

    if chromosome in chromosome_mapping:
        return chromosome_mapping[chromosome]
    else:
        return chromosome

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
    with gzip.open(mane_gff, 'rt') as fin:
        for line in fin:

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
                                mane_select_ens_id.add(trans_id)
                            if tag == "MANE_Plus_Clinical":
                                mane_clin_ens_id.add(trans_id)

    return mane_clin_ens_id, mane_select_ens_id

mane_clin_ens_id, mane_select_ens_id = get_mane(mane_gff)



def parse_gencode(gencode_gff, mane_clin_ens_id: set, mane_select_ens_id: set):
    global_gene_name_dict = {}

    current_ensembl_gene = None
    current_ensembl_transcript = None
    current_ensembl_exon = None
    current_ensembl_cds = None

    gencode_genename_genobject = dict()
    gencode_genename_geneid_transobject = dict()
    
    # Creates a dictionary with a default value of an empty dictionary
    gencode_genename_geneid_transid_exonobject = defaultdict(dict)
    gencode_genename_geneid_transid_cdsobject = defaultdict(dict)

    with gzip.open(gencode_gff, 'rt') as fin:
        feature_fields = set()
        for line in fin:
            gene_dict = {}
            line = line.rstrip("\n")

            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            feature = tmp[2].lower()
            chr = tmp[0].replace("chr", "")
            pos = tmp[3]
            end = tmp[4]
            info = tmp[8].split(";")

            # parsing genes and creating a EnsemblGene dataclass with its information
            if feature == "gene":
                if current_ensembl_gene is not None:
                    if gene_name not in gencode_genename_genobject:
                        gencode_genename_genobject[gene_name] = list()
                        gencode_genename_genobject[gene_name].append(current_ensembl_gene)
                    else:
                        gencode_genename_genobject[gene_name].append(current_ensembl_gene)

                gene_dict = dict()
                gene_dict["chromosome"] = chr
                gene_dict["start"] = pos
                gene_dict["end"] = end
                coordinates = f"{chr}:{pos}:{end}"

                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    feature_fields.add(field)

                    if field == "gene_name":
                        gene_name = result

                    elif field == "gene_id":
                        gene_id = result

                    gene_dict[field] = result
                gene_items = [
                    "gene_id",
                    "start",
                    "end",
                    "chromosome",
                    "id",
                    "tag",
                    "num_gencode_transcripts",
                    "num_refseq_transcripts",
                    "refseq_same_gene_coords",
                    "gene_name",
                    "gene_type",
                    "level",
                    "havana_gene",
                    "hgnc_id"
                ]
                filtered_gene_dict = {key : value for key, value in gene_dict.items() if key in gene_items}
                global_gene_name_dict.setdefault(gene_name, {}).setdefault("gencode", {}).setdefault(gene_id, {})["coords"] = coordinates
                current_ensembl_gene = EnsemblGenes(**filtered_gene_dict)
            
            # parsing transcripts and creating a EnsemblTranscript dataclass with its information
            elif feature == "transcript":
                if current_ensembl_transcript is not None:
                    if gene_name not in gencode_genename_geneid_transobject:
                        gencode_genename_geneid_transobject.setdefault(gene_name, dict())
                    if gene_id not in  gencode_genename_geneid_transobject:
                        gencode_genename_geneid_transobject[gene_name].setdefault(gene_id, list())
                        gencode_genename_geneid_transobject[gene_name][gene_id].append(current_ensembl_transcript)
                    else:
                        gencode_genename_geneid_transobject[gene_name][gene_id].append(current_ensembl_transcript)

                # Incrementing the attribute num_gencode_transcripts by 1 (numb of transcripts that the current gene have + 1)
                current_ensembl_gene.increment_num_gencode_transcripts()
                transcript_dict = dict()
                transcript_dict["chromosome"] = chr
                transcript_dict["start"] = pos
                transcript_dict["end"] = end
                coordinates = f"{chr}:{pos}:{end}"

                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    feature_fields.add(field)

                    if field == "transcript_id":
                        transcript_id = result

                    elif field == "gene_name":
                        gene_name = result
                    
                    elif field == "gene_id":
                        gene_id = result
                    
                    # adding all transcript information in transcript_dict
                    transcript_dict[field] = result
                    # adding mane_select and clin
                if transcript_id in mane_select_ens_id:
                    transcript_dict["mane_select"] = True
                if transcript_id in mane_clin_ens_id:
                    transcript_dict["mane_clin"] = True

                transcript_items = [
                    "transcript_id",
                    "mane_select",
                    "mane_clin",
                    "start",
                    "end",
                    "chromosome",
                    "transcript_type",
                    "transcript_name",
                    "id",
                    "havana_transcript",
                    "protein_id",
                    "gene_id",
                    "lrg_transcript",
                    "lrg_id",
                    "ccds"
                ]
                filtered_transcript_dict = {key : value for key, value in transcript_dict.items() if key in transcript_items}
                current_ensembl_transcript = EnsemblTranscript(**filtered_transcript_dict)

            # parsing exons and creating a EnsemblExon dataclass with its information
            elif feature == "exon":
                current_ensembl_transcript.increment_numb_exons()
                if current_ensembl_exon is not None:
                    if gene_name not in gencode_genename_geneid_transid_exonobject:
                        gencode_genename_geneid_transid_exonobject.setdefault(gene_name, dict())
                    if gene_id not in  gencode_genename_geneid_transid_exonobject[gene_name]:
                        gencode_genename_geneid_transid_exonobject[gene_name].setdefault(gene_id, dict())
                    if transcript_id not in gencode_genename_geneid_transid_exonobject[gene_name][gene_id]:
                        gencode_genename_geneid_transid_exonobject[gene_name][gene_id].setdefault(transcript_id, list())
                        gencode_genename_geneid_transid_exonobject[gene_name][gene_id][transcript_id].append(current_ensembl_exon)
                    else:
                        gencode_genename_geneid_transid_exonobject[gene_name][gene_id][transcript_id].append(current_ensembl_exon)
                exon_dict = dict()
                exon_dict["chromosome"] = chr
                exon_dict["start"] = pos
                exon_dict["end"] = end
                coordinates = f"{chr}:{pos}:{end}"

                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    feature_fields.add(field)

                    if field == "exon_id":
                        exon_id = result

                    elif field == "gene_name":
                        gene_name = result
                    
                    elif field == "gene_id":
                        gene_id = result

                    exon_dict[field] = result
                
                exon_items = [
                    "exon_id",
                    "id",
                    "exon_number",
                    "start",
                    "end",
                    "chromosome",
                    "transcript_id"
                ]

                filtered_exon_dict = {key : value for key, value in exon_dict.items() if key in exon_items}
                current_ensembl_exon = EnsemblExons(**filtered_exon_dict)

            # parsing genes and creating a EnsemblCds dataclass with its information
            elif feature == "cds":
                current_ensembl_transcript.increment_numb_cds()

                if current_ensembl_cds is not None:
                    if gene_name not in gencode_genename_geneid_transid_cdsobject:
                        gencode_genename_geneid_transid_cdsobject.setdefault(gene_name, dict())
                    if gene_id not in  gencode_genename_geneid_transid_cdsobject[gene_name]:
                        gencode_genename_geneid_transid_cdsobject[gene_name].setdefault(gene_id, dict())
                    if transcript_id not in gencode_genename_geneid_transid_cdsobject[gene_name][gene_id]:
                        gencode_genename_geneid_transid_cdsobject[gene_name][gene_id].setdefault(transcript_id, list())
                        gencode_genename_geneid_transid_cdsobject[gene_name][gene_id][transcript_id].append(current_ensembl_cds)
                    else:
                        gencode_genename_geneid_transid_cdsobject[gene_name][gene_id][transcript_id].append(current_ensembl_cds)

                cds_dict = dict()
                cds_dict["chromosome"] = chr
                cds_dict["start"] = pos
                cds_dict["end"] = end
                coordinates = f"{chr}:{pos}:{end}"

                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    feature_fields.add(field)

                    if field == "id":
                        cds_id = result

                    elif field == "gene_name":
                        gene_name = result
                    
                    elif field == "parent":
                        transcript_id = result
                    
                    elif field == "gene_id":
                        gene_id = result

                    cds_dict[field] = result
                cds_items = [
                    "id",
                    "start",
                    "end",
                    "chromosome",
                    "exon_id",
                    "exon_number",
                    "transcript_id"
                ]

                filtered_cds_dict = {key : value for key, value in cds_dict.items() if key in cds_items}
                current_ensembl_cds = EnsemblCds(**filtered_cds_dict)
                
    
    # Adding the last gene, transcript, exon and cds found in the gff file
    # as it is only added when the following line with the same element is found
    # as consequence the lasts components are not included in the dictionary until now
    if gene_name not in gencode_genename_geneid_transid_cdsobject:
        gencode_genename_geneid_transid_cdsobject.setdefault(gene_name, dict())
    if gene_id not in  gencode_genename_geneid_transid_cdsobject[gene_name]:
        gencode_genename_geneid_transid_cdsobject[gene_name].setdefault(gene_id, dict())
    if transcript_id not in gencode_genename_geneid_transid_cdsobject[gene_name][gene_id]:
        gencode_genename_geneid_transid_cdsobject[gene_name][gene_id].setdefault(transcript_id, list())
        gencode_genename_geneid_transid_cdsobject[gene_name][gene_id][transcript_id].append(current_ensembl_cds)

    else:
        gencode_genename_geneid_transid_cdsobject[gene_name][gene_id][transcript_id].append(current_ensembl_cds)

    if gene_name not in gencode_genename_geneid_transid_exonobject:
        gencode_genename_geneid_transid_exonobject.setdefault(gene_name, dict())
    if gene_id not in  gencode_genename_geneid_transid_exonobject[gene_name]:
        gencode_genename_geneid_transid_exonobject[gene_name].setdefault(gene_id, dict())
    if transcript_id not in gencode_genename_geneid_transid_exonobject[gene_name][gene_id]:
        gencode_genename_geneid_transid_exonobject[gene_name][gene_id].setdefault(transcript_id, list())
        gencode_genename_geneid_transid_exonobject[gene_name][gene_id][transcript_id].append(current_ensembl_exon)
    else:
        gencode_genename_geneid_transid_exonobject[gene_name][gene_id][transcript_id].append(current_ensembl_exon)


    if gene_name not in gencode_genename_geneid_transobject:
        gencode_genename_geneid_transobject.setdefault(gene_name, dict())
    if gene_id not in  gencode_genename_geneid_transobject:
        gencode_genename_geneid_transobject[gene_name].setdefault(gene_id, list())
        gencode_genename_geneid_transobject[gene_name][gene_id].append(current_ensembl_transcript)
    else:
        gencode_genename_geneid_transobject[gene_name][gene_id].append(current_ensembl_transcript)

    gencode_genename_genobject[gene_name] = current_ensembl_gene
        
    return (gencode_genename_genobject, gencode_genename_geneid_transobject, gencode_genename_geneid_transid_exonobject, gencode_genename_geneid_transid_cdsobject)


(
    gencode_genename_genobject,
    gencode_genename_geneid_transobject,
    gencode_genename_geneid_transid_exonobject,
    gencode_genename_geneid_transid_cdsobject
) = parse_gencode(gencode_gff, mane_clin_ens_id, mane_select_ens_id)

print(gencode_genename_geneid_transid_cdsobject)