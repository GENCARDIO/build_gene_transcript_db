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

# class Gene:
#     def __init__(self, gene_name, gene_coords, refseq_id = None, gencode_id = None):
#         self.gene_name = gene_name
#         self.refseq_id = refseq_id
#         self.gencode_id = gencode_id
#         self.gene_coords = gene_coords

# class Transcript(Gene):
#     def __init__(self, transcript_coords, gencode_transcript_id = None, refseq_transcript_id = None):
#         self.transcript_coords = transcript_coords
#         self.gencode_transcript_id = 


def parse_gencode(gencode_gff, global_gene_name_dict):

    # Creates a dictionary with a default value of an empty dictionary
    gencode_gene_dict = defaultdict(dict)
    gencode_transcript_dict = defaultdict(dict)
    gencode_exon_dict = defaultdict(dict)
    gencode_cds_dict = defaultdict(dict)

    with gzip.open(gencode_gff, 'rt') as fin:
        feature_fields = set()
        for line in fin:
            line = line.rstrip("\n")

            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            feature = tmp[2].lower()
            chr = tmp[0].replace("chr", "")
            pos = tmp[3]
            end = tmp[4]
            info = tmp[8].split(";")

            if feature == "gene":
                gene_dict = dict()
                gene_dict["chromosome"] = chr
                gene_dict["start"] = pos
                gene_dict["end"] = end
                coordinates = f"{chr}:{pos}:{end}"

                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    feature_fields.add(field)

                    # if field == "gene_name":
                    #     if not result in gene_name:
                    #         gene_name.append(result)
                    #     else:
                    #         same_gene_name += 1

                    if field == "gene_name":
                        gene_name = result

                    elif field == "gene_id":
                        gene_id = result

                    # if field == "ID":
                        # there are genes with same identifier but with _PAR_Y
                        # at the final of the ID e.g. ENSG00000229232.6_PAR_Y
                        # (that we will differenciate between them as they are
                        # pseudoautosomal region on the Y chromosome
                        # if "PAR_Y" not in result.split(".")[1]:
                        #     # not taking into account gene version
                        #     result = result.split(".")[0]
                        # else:
                        #     result = str(result.split(".")[0]) + "_PAR_Y"
                    gene_dict[field] = result
                global_gene_name_dict.setdefault(gene_name, {}).setdefault("gencode", {}).setdefault(gene_id, {})["coords"] = coordinates
                gencode_gene_dict[gene_id] = gene_dict

            elif feature == "transcript":
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

                    transcript_dict[field] = result
                
                global_gene_name_dict.setdefault(gene_name, {}).setdefault("gencode", {}).setdefault(gene_id, {}).setdefault("transcripts", {}).setdefault(transcript_id, {})["coords"] = coordinates
                gencode_transcript_dict[transcript_id] = transcript_dict

            elif feature == "exon":
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
                
                global_gene_name_dict.setdefault(gene_name, {}).setdefault("gencode", {}).setdefault(gene_id, {}).setdefault("transcripts", {}).setdefault(transcript_id, {}).setdefault("exons", {})[exon_id] = coordinates
                gencode_exon_dict[exon_id] = exon_dict

            elif feature == "cds":
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
                
                global_gene_name_dict.setdefault(gene_name, {}).setdefault("gencode", {}).setdefault(gene_id, {}).setdefault("transcripts", {}).setdefault(transcript_id, {}).setdefault("cds", {})
                if cds_id not in global_gene_name_dict[gene_name]["gencode"][gene_id]["transcripts"][transcript_id]["cds"]:
                    global_gene_name_dict.setdefault(gene_name, {}).setdefault("gencode", {}).setdefault(gene_id, {}).setdefault("transcripts", {}).setdefault(transcript_id, {}).setdefault("cds", {})[cds_id] = [coordinates]
                else:
                    global_gene_name_dict[gene_name]["gencode"][gene_id]["transcripts"][transcript_id]["cds"][cds_id].append(coordinates)
                    
                gencode_cds_dict[cds_id] = cds_dict
        
    return (global_gene_name_dict, gencode_gene_dict, gencode_transcript_dict, gencode_exon_dict, gencode_cds_dict)

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

refseq_gff = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/refseq/GRCh38_latest_genomic.gff.gz"
def get_refseq_transcripts(refseq_gff: str, global_gene_name_dict):
    '''
        Function that reads a RefSeq gff3 file and returns two structured dictionaries
            with genes and transcripts
        :param str mane_gff: RefSeq gff3 file
        :return: two dicts, one with primary keys pointing to an enst_id, and a second one
            with gene id's as primary keys
        :rtype: Two defaultdict(dict)
    '''
    gene_numb = 0
    hgnc_numb = 0
    refseq_transcripts = defaultdict(dict)
    refseq_genes_dict  = defaultdict(dict)
    refseq_cds_dict = defaultdict(dict)
    feature_dict = dict()
    feature_set = set()
    transcripts_ids = set()

    with gzip.open(refseq_gff, 'rt') as fin:
        for line in fin:
            fields = []
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
            feature_dict["chromosome"] = chr
            feature_dict["start"] = pos
            feature_dict["end"] = end
            coordinates = f"{chr}:{pos}:{end}"

            genes_types = {'match', 'repeat_region', 'tandem_repeat', 'origin_of_replication', 'transcriptional_cis_regulatory_region', 'recombination_feature', 'sequence_comparison', 'imprinting_control_region', 'chromosome_breakpoint', 'dnasei_hypersensitive_site', 'insulator', 'enhancer', 'caat_signal', 'epigenetically_modified_region', 'd_loop', 'protein_binding_site', 'minisatellite', 'meiotic_recombination_region', 'nucleotide_cleavage_site', 'region', 'replication_start_site', 'regulatory_region', 'nucleotide_motif', 'tata_box', 'conserved_region', 'non_allelic_homologous_recombination_region', 'sequence_alteration_artifact', 'microsatellite', 'gc_rich_promoter_region', 'gene', 'response_element', 'cdna_match', 'centromere', 'silencer', 'mobile_genetic_element', 'mitotic_recombination_region', 'sequence_secondary_structure', 'matrix_attachment_site', 'replication_regulatory_region', 'sequence_alteration', 'dispersed_repeat', 'locus_control_region', 'cage_cluster', 'direct_repeat', 'repeat_instability_region', 'enhancer_blocking_element', 'sequence_feature', 'promoter', 'pseudogene'}
            if feature in genes_types:
                gene_id = None
                
                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    if field == "dbxref":
                        dbxref_results = result.split(",")
                        for dbxref_result in dbxref_results:
                            if "GeneID" in dbxref_result:
                                gene_id = dbxref_result.split(":")[1]
                                feature_dict["id"] = gene_id   

                    elif field == "name":
                        gene_name = result

                global_gene_name_dict.setdefault(gene_name, {}).setdefault("refseq", {}).setdefault(gene_id, {})[coordinates] = {}

            transcripts_types = ["transcript", "primary_transcript", "lnc_rna", "mrna", "snrna"]
            if feature in transcripts_types:
                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    if field == "parent":
                        if "-" in result:
                        # removing gene- from the gene name
                            gene_name = result.split("-", 1)[-1]
                    elif field == "id":
                        transcript_id = result
                        transcripts_ids.add(transcript_id)
                    elif field == "dbxref":
                        dbxref_results = result.split(",")
                        for dbxref_result in dbxref_results:
                            if "GeneID" in dbxref_result:
                                gene_id = dbxref_result.split(":")[1]
                                feature_dict["id"] = gene_id

                            # elif "hgnc" in dbxref_result:
                            #     hgnc_id = dbxref_result.split(":")[-1]
                            #     hgnc_id = f"HGNC:{hgnc_id}"
                            #     feature_dict["hgnc_id"] = hgnc_id
                            #     hgnc_numb += 1
                global_gene_name_dict.setdefault(gene_name, {}).setdefault("refseq", {}).setdefault(gene_id, {}).setdefault("transcripts", {}).setdefault(transcript_id, {})["coords"] = coordinates

            if feature == "exon":

                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    if field == "id":
                        exon_id = result
                    elif field == "gene":
                        gene_name = result
                    elif field == "parent":
                        transcript_id = result
                    elif field == "dbxref":
                        dbxref_results = result.split(",")
                        for dbxref_result in dbxref_results:
                            if "GeneID" in dbxref_result:
                                gene_id = dbxref_result.split(":")[1]
                                feature_dict["id"] = gene_id
                            # elif "mirbase" in dbxref_result:
                            #     mirbase_id = dbxref_result.split(":")[1]
                            #     feature_dict["mirbase_id"] = mirbase_id
                            # elif "mim" in dbxref_result:
                            #     mim_id = dbxref_result.split(":")[1]
                            #     feature_dict["mim_id"] = mim_id
                    elif field == "gene_biotype":
                        field = "gene_type"
                    elif field == "gene":
                        id = result
                    feature_dict[field] = result
                if not transcript_id.startswith("rna-MIR"):
                    global_gene_name_dict.setdefault(gene_name, {}).setdefault("refseq", {}).setdefault(gene_id, {}).setdefault("transcripts", {}).setdefault(transcript_id, {}).setdefault("exons", {})[exon_id] = coordinates

                # refseq_genes_dict[id] = feature_dict
            
            
            if feature == "cds":

                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    if field == "id":
                        cds_id = result
                    elif field == "gene":
                        gene_name = result
                    elif field == "parent":
                        transcript_id = result
                    elif field == "dbxref":
                        dbxref_results = result.split(",")
                        for dbxref_result in dbxref_results:
                            if "GeneID" in dbxref_result:
                                gene_id = dbxref_result.split(":")[1]
                                feature_dict["id"] = gene_id
                global_gene_name_dict.setdefault(gene_name, {}).setdefault("refseq", {}).setdefault(gene_id, {}).setdefault("transcripts", {}).setdefault(transcript_id, {}).setdefault("cds", {})
                if cds_id not in global_gene_name_dict[gene_name]["refseq"][gene_id]["transcripts"][transcript_id]["cds"]:
                    global_gene_name_dict.setdefault(gene_name, {}).setdefault("refseq", {}).setdefault(gene_id, {}).setdefault("transcripts", {}).setdefault(transcript_id, {}).setdefault("cds", {})[cds_id] = [coordinates]
                else:
                    global_gene_name_dict[gene_name]["refseq"][gene_id]["transcripts"][transcript_id]["cds"][cds_id].append(coordinates)
                    
    return refseq_genes_dict, global_gene_name_dict, transcripts_ids

global_gene_name_dict = defaultdict(dict)
refseq_genes_dict, global_gene_name_dict, transcripts_ids = get_refseq_transcripts(refseq_gff, global_gene_name_dict)

(   
    global_gene_name_dict,
    gencode_gene_dict,
    gencode_transcript_dict,
    gencode_exon_dict,
    gencode_cds_dict
) = parse_gencode(gencode_gff, global_gene_name_dict)

json_dict_file = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/gene_dict.json"
with open(json_dict_file, "w") as file:
    json.dump(global_gene_name_dict, file)






mane_gff = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/mane/MANE.GRCh38.v1.1.ensembl_genomic.gff.gz"
def get_mane(mane_gff: str, gencode_transcript_dict, gencode_exon_dict):
    '''
        Function that reads a Mane gff3 file and returns two structured dictionaries
            with genes and transcripts
        :param str mane_gff: RefSeq gff3 file
        :return: two dicts, one with primary keys pointing to an enst_id, and a second one
            with gene id's as primary keys
        :rtype: Two defaultdict(dict)
    '''
    features = set()
    
    with gzip.open(mane_gff, 'rt') as fin:
        for line in fin:
            mane_select = False
            mane_clin = False
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
            coordinates = f"{chr}:{pos}:{end}"
            
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
                if trans_id in gencode_transcript_dict:
                    if mane_select is True:
                        gencode_transcript_dict[trans_id]["mane_select"] = True
                    else:
                        gencode_transcript_dict[trans_id]["mane_select"] = None
                    if mane_clin is True:
                        gencode_transcript_dict[trans_id]["mane_clin"] = True
                    else:
                        gencode_transcript_dict[trans_id]["mane_clin"] = None
                
            if feature == "exon":
                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    if field == "exon_id":
                        exon_id = result
                    if field == "tag":
                        tags = result.split(",")
                        for tag in tags:
                            features.add(tag)
                            if tag == "MANE_Select":
                                mane_select = True
                            if tag == "MANE_Plus_Clinical":
                                mane_clin = True
                if exon_id in gencode_exon_dict:
                    if mane_select is True:
                        gencode_exon_dict[exon_id]["mane_select"] = True
                    else:
                        gencode_exon_dict[exon_id]["mane_select"] = False
                    if mane_clin is True:
                        gencode_exon_dict[exon_id]["mane_clin"] = True
                    else:
                        gencode_exon_dict[exon_id]["mane_clin"] = False

    return(gencode_transcript_dict, gencode_exon_dict)

gencode_transcript_dict, gencode_exon_dict = get_mane(mane_gff, gencode_transcript_dict, gencode_exon_dict)

def get_lrg_trancripts(lrg_txt: str, gencode_transcript_dict: dict):

    with open(lrg_txt) as f:
        for line in f:
            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            lrg_id = tmp[0]
            hgnc_symbol = tmp[1]
            refseq_genomic_id = tmp[2]
            lrg_transcript = tmp[3]
            refseq_transcript_id = tmp[4]
            ensembl_trans_id = tmp[5]
            ccds = tmp[6]

            if ensembl_trans_id in gencode_transcript_dict:
                gencode_transcript_dict[ensembl_trans_id]["lrg_id"] = lrg_id
                gencode_transcript_dict[ensembl_trans_id]["lrg_transcript"] = lrg_transcript
                gencode_transcript_dict[ensembl_trans_id]["ccds"] = ccds

    return gencode_transcript_dict

gencode_transcript_dict = get_lrg_trancripts("/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/lrg/list_LRGs_transcripts_xrefs.txt", gencode_transcript_dict)

def compare_gene_composition():
    """
    Global_gene_name_dict has the following structure:
    global_gene_name_dict[gene_name]
    Then it seperats if the genes provide from refseq or gencode
    global_gene_name_dict[gene_name]["refseq"]
    global_gene_name_dict[gene_name]["gencode"]
    Then it seperates between genes and transcripts (just following the structure of gencode but refseq has the same):
    global_gene_name_dict[gene_name]["gencode"][gene_id]["coords"] = coords
    global_gene_name_dict[gene_name]["gencode"]["transcripts"][transcript_id]["coords"] = coords
    Finally for each transcript we have its exons and cds coordinates
    global_gene_name_dict[gene_name]["gencode"]["transcripts"][transcript_id]["exons"][exon_id] = coords_list
    global_gene_name_dict[gene_name]["gencode"]["transcripts"][transcript_id]["cds"][cds_id] = coords_list

    From here we iterate on over gene_name to detect the matches between refseq and gencode equivalent:
        genes coords
        transcripts coords
        exon coords that are composing each transcript
        cds coords that are composing each transcript
    
    Return:
        geneid_stats: has the following structure:
        geneid_stats[gencode_id]["equal_gene_coords"] = gene_refseq_id
        geneid_stats[gencode_transcript_id]["equal_trans_coords"] = transcript_refseq_id
        geneid_stats[gencode_transcript_id]["equal_exons"] = transcript_refseq_id
        geneid_stats[gencode_transcript_id]["equal_cds"] = transcript_refseq_id

    """
    with open(json_dict_file, "r") as file:
        global_gene_name_dict = json.load(file)
    geneid_stats = {}

    for gene_name, gene_info in global_gene_name_dict.items():
        refseq_coords = list()
        refseq_gene_coords = dict()
        # If gene_name has associated a gencode gene and refseq gene
        if "gencode" in gene_info and "refseq" in gene_info:
            gencode_gene_info = gene_info["gencode"]
            refseq_gene_info = gene_info["refseq"]
            
            for refseq_gene_id, refseq_gene_item in refseq_gene_info.items():
                # Obtaining all refseq coordinates for a given gene_name
                for coord, _ in refseq_gene_item.items():
                    refseq_coords.append(coord)
                    refseq_gene_coords[coord] = refseq_gene_id
            # Obtaining all gencode coordinates for a given gene_name
            for gencode_gene_id, gencode_gene_item in gencode_gene_info.items():
                geneid_stats.setdefault("gene", {})
                geneid_stats["gene"][gencode_gene_id] = {}
                coord = gencode_gene_item["coords"]
                # comparing if gencode coord of the gene are the same of the refseq coords for the same gene name
                if coord in refseq_coords:
                    geneid_stats["gene"][gencode_gene_id]["refseq_same_gene_coords"] = refseq_gene_coords[coord]
                else:
                    geneid_stats["gene"][gencode_gene_id]["refseq_same_gene_coords"] = None

                # give transcripts info into gencode_transcripts and refseq_transcripts variables:
                if "transcripts" in gencode_gene_item:
                    gencode_transcripts = gencode_gene_item["transcripts"]
                else:
                    gencode_transcripts = None
                if "transcripts" in refseq_gene_item:
                    refseq_transcripts = refseq_gene_item["transcripts"]
                else:
                    refseq_transcripts = None
                
                # Creating dictionaries where we will store all the transcripts, exons and cds
                # refseq coordinates for a given gene name to be compared with gencode coords
                refseq_trans_coords = dict()
                refseq_exons_coords = dict()
                refseq_cds_coords = dict()
                
                geneid_stats["gene"][gencode_gene_id]["num_gencode_transcripts"] = len(gencode_transcripts)

                #if refseq and gencode have transcripts associated for a given gene_name
                if refseq_transcripts and gencode_transcripts:
                    geneid_stats["gene"][gencode_gene_id]["num_refseq_transcripts"] = len(refseq_transcripts)

                    # Obtaining refseq exons and cds coords for a given gene_name
                    for refseq_transcript_id, refseq_transcript_info in refseq_transcripts.items():
                        if "coords" in refseq_transcript_info:
                            coords = refseq_transcript_info["coords"]
                            if "cds" in refseq_transcript_info and "exons" not in refseq_transcript_info:
                                print(transcript_id, "has cds but not exons!")
                            if "exons" in refseq_transcript_info:
                                exons_info = refseq_transcript_info["exons"]
                                for exon_id, exon_coord in exons_info.items():
                                    if refseq_transcript_id not in refseq_exons_coords:
                                        #In this variable will be stored a set of the exon coordinates of the transcript
                                        refseq_exons_coords[refseq_transcript_id] = {exon_coord}
                                    else:
                                        refseq_exons_coords[refseq_transcript_id].add(exon_coord)
                                refseq_exons_coords[refseq_transcript_id]
                            if "cds" in refseq_transcript_info:
                                cds_info = refseq_transcript_info["cds"]
                                for cds_id, cds_coords_list in cds_info.items():
                                    if refseq_transcript_id not in refseq_cds_coords:
                                        # In this variable is stored a set of cds coords of the transcript
                                       refseq_cds_coords[refseq_transcript_id] = set(cds_coords_list)
                            if coords not in refseq_trans_coords:
                                refseq_trans_coords[coords] = [refseq_transcript_id]
                            else:
                                refseq_trans_coords[coords].append(refseq_transcript_id)

                    # Obtaining gencode exons and cds coords for a given gene_name and compare the coordinates with refseq
                    for transcript_id, gencode_transcript_info in gencode_transcripts.items():
                        gencode_exon_set = set()
                        gencode_cds_set = set()
                        if "coords" in gencode_transcript_info:
                            trans_coord = gencode_transcript_info["coords"]
                            if trans_coord in refseq_trans_coords:
                                geneid_stats.setdefault("transcript", {}).setdefault(transcript_id, {})["refseq_same_trans_coords"] = refseq_trans_coords[trans_coord]
                            else:
                                geneid_stats.setdefault("transcript", {}).setdefault(transcript_id, {})["refseq_same_trans_coords"] = []
                        if "exons" in gencode_transcript_info:
                            exons_info = gencode_transcript_info["exons"]
                            for gencode_exon_id, gencode_exon_coord in exons_info.items():
                                gencode_exon_set.add(gencode_exon_coord)
                            for refseq_transcript_id, refseq_exon_set in refseq_exons_coords.items():
                                # comparing each gencode trancript exon composition with
                                # every refseq transcript exon composition
                                if gencode_exon_set == refseq_exon_set:
                                    if transcript_id not in geneid_stats["transcript"]:
                                        geneid_stats.setdefault("transcript", {}).setdefault(transcript_id, {})
                                    if "refseq_same_exons" not in geneid_stats["transcript"][transcript_id]:
                                        geneid_stats["transcript"][transcript_id].setdefault("refseq_same_exons", {})
                                    geneid_stats["transcript"][transcript_id]["refseq_same_exons"][refseq_transcript_id] = True

                                else:
                                    geneid_stats.setdefault("transcript").setdefault(transcript_id, {})
                        if "cds" in gencode_transcript_info:
                            cds_info = gencode_transcript_info["cds"]
                            for gencode_cds_id, gencode_cds_coord in cds_info.items():
                                gencode_cds_set = set(gencode_cds_coord)
                            for refseq_transcript_id, refseq_cds_set in refseq_cds_coords.items():
                                # comparing each gencode transcript cds composition with
                                # every refseq transcript cds composition
                                if gencode_cds_set == refseq_cds_set:
                                    print(gencode_cds_set, refseq_cds_set)
                                    if "refseq_same_cds" not in geneid_stats["transcript"][transcript_id]:
                                        geneid_stats["transcript"][transcript_id].setdefault("refseq_same_cds", {})
                                    geneid_stats["transcript"][transcript_id]["refseq_same_cds"][refseq_transcript_id] = True


    return (geneid_stats)



geneid_stats = compare_gene_composition()


def add_geneid_gene_stats_to_global_dict(geneid_stats, gencode_gene_dict):

    for gencode_gene_id, items in geneid_stats["gene"].items():
        if gencode_gene_dict[gencode_gene_id]:
            gencode_gene_dict[gencode_gene_id].update(items)
    
    return gencode_gene_dict

def add_geneid_trans_stats_to_global_dict(geneid_stats, gencode_transcript_dict):
    for gencode_transcript_id, items in geneid_stats["transcript"].items():
        if gencode_transcript_dict[gencode_transcript_id]:
            print(items, "item")
            gencode_transcript_dict[gencode_transcript_id].update(items)
    
    return gencode_transcript_dict





gencode_gene_dict = add_geneid_gene_stats_to_global_dict(geneid_stats, gencode_gene_dict)
gencode_transcript_dict = add_geneid_trans_stats_to_global_dict(geneid_stats, gencode_transcript_dict)

def create_database():
    mysql_psw = os.environ["MYSQL_PWD"]
    # Now define and create the genes database
    genes_db_name = "genes_db"
    engine  = create_engine(f'mysql+mysqlconnector://ocanal:{mysql_psw}@localhost:3306')
    Session = sessionmaker(bind=engine)
    session = Session()
    Base = declarative_base()

    try:
        engine.execute(f"CREATE DATABASE IF NOT EXISTS {genes_db_name}")
        engine.execute(f"USE {genes_db_name}")
        print("database genes_db created successfully")
    
    except exc.SQLAlchemyError as e:
        print("Error ocurred while creating the database:", str(e))

    class Genes(Base):
        __tablename__ = "Genes"
        gene_id = Column(String(50), primary_key=True)
        start = Column(Integer, nullable=False)
        end = Column(Integer, nullable=False)
        chromosome = Column(String(20), nullable=False)
        id = Column(String(50), nullable=False)
        tag = Column(String(100))
        num_gencode_transcripts = Column(Integer)
        num_refseq_transcripts = Column(Integer)
        refseq_same_gene_coords = Column(String(50))
        gene_name = Column(String(50))
        gene_type = Column(String(50))
        level = Column(Integer)
        havana_gene = Column(String(50))
        hgnc_id = Column(String(50))
        
        transcripts = relationship("Transcripts", back_populates="gene")

        def __repr__(self):
            return f"Gene {self.gene_id}"

    association_table_trans_coords = Table(
        'transcripts_refseq_trans_coords',
        Base.metadata,
        Column('transcript_id', String(50), ForeignKey('Transcripts.transcript_id')),
        Column('refseq_id', String(50), ForeignKey('RefseqTranscripts.refseq_id'))
    )

    association_table_exons = Table(
        'transcripts_refseq_exons',
        Base.metadata,
        Column('transcript_id', String(50), ForeignKey('Transcripts.transcript_id')),
        Column('refseq_id', String(50), ForeignKey('RefseqTranscripts.refseq_id'))
    )

    association_table_cds = Table(
        'transcripts_refseq_cds',
        Base.metadata,
        Column('transcript_id', String(50), ForeignKey('Transcripts.transcript_id')),
        Column('refseq_id', String(50), ForeignKey('RefseqTranscripts.refseq_id'))
    )


    class RefseqTranscripts(Base):
        __tablename__ = "RefseqTranscripts"
        refseq_id = Column(String(50), primary_key=True)

        def __repr__(self):
            return f"Refseq Transcript{self.refseq_id}"


    class Transcripts(Base):
        __tablename__ = "Transcripts"
        transcript_id = Column(String(50), primary_key=True)
        mane_select = Column(Boolean, nullable=True, default=None)
        mane_clin = Column(Boolean, nullable=True, default=None)
        start = Column(Integer, nullable=False)
        end = Column(Integer, nullable=False)
        chromosome = Column(String(10), nullable=False)
        transcript_type = Column(String(50))
        transcript_name = Column(String(50))
        id = Column(String(50), nullable=False)
        havana_transcript = Column(String(50))
        protein_id = Column(String(50))
        gene_id = Column(String(50), ForeignKey("Genes.gene_id"))
        lrg_transcript = Column(String(50))
        lrg_id = Column(String(50))
        ccds = Column(String(50))


        # Relationships
        gene = relationship("Genes", back_populates="transcripts")
        exons = relationship("Exons", back_populates="transcript")
        refseq_same_trans_coords = relationship(
            "RefseqTranscripts",
            secondary=association_table_trans_coords,
            backref="transcripts",
        )
        refseq_same_exons = relationship(
            "RefseqTranscripts",
            secondary=association_table_exons,
            backref="transcripts_exons",
        )
        refseq_same_cds = relationship(
            "RefseqTranscripts",
            secondary=association_table_cds,
            backref="transcripts_cds",
        )

        def __repr__(self):
            return f"Transcript {self.transcript_id}"


    class Exons(Base):
        __tablename__ = "Exons"
        exon_id = Column(String(50), primary_key=True)
        id = Column(String(50), nullable=False)
        exon_number = Column(Integer)
        start = Column(Integer, nullable=False)
        end = Column(Integer, nullable=False)
        chromosome = Column(String(50), nullable=False)
        transcript_id = Column(String(50), ForeignKey("Transcripts.transcript_id"))

        #Relationships
        transcript = relationship("Transcripts", back_populates="exons")
        cds = relationship("Cds", back_populates="exon")

        def __repr__(self):
            return f"Exon {self.exon_id}"

    class Cds(Base):
        __tablename__ = "Cds"
        id = Column(String(50), primary_key=True)
        start = Column(Integer, nullable=False)
        end = Column(Integer, nullable=False)
        chromosome = Column(String(50), nullable=False)
        exon_id = Column(String(50), ForeignKey("Exons.exon_id"))

        # Relationships
        exon = relationship("Exons", back_populates="cds")

        def __repr__(self):
            return f"CDS: {self.id}"
    
    # inspector = inspect(engine)
    # # Check if the Genes table exists in the database
    # if not inspector.has_table("Genes"):
    #     # Create all tables defined in the base class
    #     Base.metadata.create_all(engine)

    return (
        Genes,
        Transcripts,
        RefseqTranscripts,
        Exons,
        Cds,
        association_table_trans_coords,
        association_table_exons,
        association_table_cds,
        session
    )


(
    Genes,
    Transcripts,
    RefseqTranscripts,
    Exons,
    Cds,
    association_table_trans_coords,
    association_table_exons,
    association_table_cds,
    session
) = create_database()


def insert_db_gencode_gene_data(gencode_gene_dict, session, Genes):
    """
    From a dictionary with [gene_id]:[feature]:[values]
    inserts the features and values into the genes table

    Params:
        gencode_gene_dict[gene_id]:[feature]:[values]
        session: session of the database, returned by function create_database
    """
    valid_gene_keys = {
        'gene_id',
        'start',
        'end',
        'chromosome',
        'id',
        'tag',
        'gene_name',
        'gene_type',
        'level',
        "havana_gene",
        "hgnc_id",
        "refseq_same_gene_coords",
        "num_gencode_transcripts",
        "num_refseq_transcirpts"
    }

    for gene, gene_items in gencode_gene_dict.items():
        filtered_gene_data = {}
        
        # Filtering the dict to only contain the key-value pairs that are contained in the database
        filtered_gene_data = {key: value for key, value in gene_items.items() if key in valid_gene_keys}

        # Using ** it recieves the following data: gene_id="gene_id", start=100, end=200 ...
        new_gene = Genes(**filtered_gene_data)

        try:
            print(f"inserting gene {new_gene}")
            # Add new gene to the session
            session.add(new_gene)
            session.commit()
        except:
            msg = "Error in database writing"
            raise SQLAlchemyError(msg)
        else:
            msg = f"gene {new_gene} record successfully"
            logging.info(msg)

def insert_refseq_transcripts_data(transcripts_ids, session, RefseqTranscripts):
    for transcript_id in transcripts_ids:
        new_refseq_transcript = RefseqTranscripts(refseq_id=transcript_id)
        try:
            print(f"inserting refseq transcript {transcript_id}")
            session.add(new_refseq_transcript)
            session.commit()
        except:
            msg = "Error in refseq transcript database writing"
            raise SQLAlchemyError(msg)
        

def insert_db_gencode_transcirpts_data(gencode_transcript_dict, session, Transcripts, RefseqTranscripts, association_table_trans_coords, association_table_exons, association_table_cds):
    """
    From a dictionary with [transcript_id]:[feature]:[values]
    inserts the features and values into the transcripts table

    Params:
        gencode_transcript_dict[transcript_id]:[feature]:[values]
        session: session of the database, returned by function create_database
        Transcripts: Transcripts database table
    """

    valid_transcript_keys = {
        'transcript_id',
        'mane_select',
        'mane_clin',
        'start',
        'end',
        'chromosome',
        'transcript_type',
        'transcript_name',
        'id',
        'havana_transcript',
        'protein_id',
        "gene_id",
        "refseq_same_trans_coords",
        "refseq_same_exons",
        "refseq_same_cds",
        "lrg_transcript",
        "lrg_id",
        "ccds"
    }

    for transcript, transcript_items in gencode_transcript_dict.items():
        refseq_transcript_db_object = None
        filtered_transcript_data = {}
        db_object_list_refseq_same_trans = []
        db_object_list_refseq_same_exons = []
        db_object_list_refseq_same_cds = []
        # Assigning db refseq objects to the dict items:"refseq_same_trans_coords",
        # "refseq_same_exons","refseq_same_cds"
        if "refseq_same_trans_coords" in transcript_items:
            if transcript_items["refseq_same_trans_coords"] != []:
                for transcript_id in transcript_items["refseq_same_trans_coords"]:
                    # print(f"refseqid: {transcript_id}")
                    # # checking if the refseq_id is already in the database:
                    # exists_query = session.query(exists().where(RefseqTranscripts.refseq_id == transcript_id))
                    # refseq_exists = session.query(exists_query.subquery()).scalar()
                    # # if refseq_id does not exists in the database:
                    # if not refseq_exists:
                    #     print(transcript_id, "refseq_id should only contain 1 item")
                    #     print(len(transcript_id))
                    #     refseq_db_object = RefseqTranscripts(refseq_id=transcript_id)
                    #     print(refseq_db_object)
                    # else:
                    refseq_transcript_db_object = session.query(RefseqTranscripts).filter_by(refseq_id=transcript_id).first()
                    db_object_list_refseq_same_trans.append(refseq_transcript_db_object)
                transcript_items["refseq_same_trans_coords"] = db_object_list_refseq_same_trans
        #         transcript_items["refseq_same_trans_coords"] = db_object_list_refseq_same_trans
        #     else:
        #         transcript_items["refseq_same_trans_coords"] = []
        # else:
        #     transcript_items["refseq_same_trans_coords"] = []
        
        if "refseq_same_exons" in transcript_items:
            for transcript_id, _ in transcript_items["refseq_same_exons"].items():
                # checking if the refseq_id is already in the database:
                exists_query = session.query(exists().where(RefseqTranscripts.refseq_id == transcript_id))
                refseq_exists = session.query(exists_query.subquery()).scalar()
                # if refseq_id does not exists in the database:
                # if not refseq_exists:
                #     refseq_db_exon_object = RefseqTranscripts(refseq_id=transcript_id)
                # else:
                refseq_db_exon_object = session.query(RefseqTranscripts).filter_by(refseq_id=transcript_id).first()

                db_object_list_refseq_same_exons.append(refseq_db_exon_object)
            transcript_items["refseq_same_exons"] = db_object_list_refseq_same_exons 
        #     transcript_items["refseq_same_exons"] = db_object_list_refseq_same_exons
        # else:
        #     transcript_items["refseq_same_exons"] = []

        if "refseq_same_cds" in transcript_items:
            for transcript_id in transcript_items["refseq_same_cds"]:
                # # checking if the refseq_id is already in the database:
                # exists_query = session.query(exists().where(RefseqTranscripts.refseq_id == transcript_id))
                # refseq_exists = session.query(exists_query.subquery()).scalar()
                # if refseq_id does not exists in the database:
                # if not refseq_exists:
                #     refseq_db_cds_object = RefseqTranscripts(refseq_id=transcript_id)
                # else:
                refseq_db_cds_object = session.query(RefseqTranscripts).filter_by(refseq_id=transcript_id).first()

                db_object_list_refseq_same_cds.append(refseq_db_cds_object)
            transcript_items["refseq_same_cds"] = db_object_list_refseq_same_cds
        #     transcript_items["refseq_same_cds"] = db_object_list_refseq_same_cds
        # else:
        #     transcript_items["refseq_same_cds"] = []


        # Filtering the dict to only contain the key-value pairs that are contained in the database
        filtered_transcript_data = {key: value for key, value in transcript_items.items() if key in valid_transcript_keys}
        # # Using ** it recieves the following data: gene_id="gene_id", start=100, end=200 ...
        new_transcript = Transcripts(**filtered_transcript_data)

        try:
            print(f"inserting transcript {new_transcript}")
            # Add new gene to the session
            session.add(new_transcript)
            session.commit()
        except:
            msg = "Error in database writing"
            raise SQLAlchemyError(msg)
        else:
            msg = f"INFO: gene {new_transcript} record successfully"
            logging.info(msg)
        
        for refseq_same_trans_coords in db_object_list_refseq_same_trans:
            association = association_table_trans_coords.insert().values(transcript_id=new_transcript.transcript_id, refseq_id=refseq_same_trans_coords.refseq_id)
            session.execute(association)
        try:
            session.commit()
        except:
            msg = "Error commiting into the database"
            logging.critical(msg)

def insert_db_gencode_exons_data(gencode_exon_dict, session, Exons):
    """
    From a dictionary with [exon_id]:[feature]:[values]
    inserts the features and values into the Exons table

    Params:
        gencode_exon_dict[exon_id]:[feature]:[values]
        session: session of the database, returned by function create_database
        Exons: Exons database table
    """


    valid_exon_keys = {
        'exon_id',
        'start',
        'end',
        'chromosome',
        'exon_number',
        'id',
        'transcript_id'
    }

    for expon, exon_items in gencode_exon_dict.items():
        filtered_exon_data = {}
        
        # Filtering the dict to only contain the key-value pairs that are contained in the database
        filtered_exon_data = {key: value for key, value in exon_items.items() if key in valid_exon_keys}

        # Using ** it recieves the following data: gene_id="gene_id", start=100, end=200 ...
        new_exon = Exons(**filtered_exon_data)

        try:
            print(f"inserting exon {new_exon}")
            # Add new gene to the session
            session.add(new_exon)
            session.commit()
        except:
            msg = "Error in database writing"
            raise SQLAlchemyError(msg)
        else:
            msg = f"INFO: gene {new_exon} record successfully"
            logging.info(msg)


def insert_db_gencode_cds_data(gencode_cds_dict, session, Cds):
    """
    From a dictionary with [cds_id]:[feature]:[values]
    inserts the features and values into the Exons table

    Params:
        gencode_cds_dict[cds_id]:[feature]:[values]
        session: session of the database, returned by function create_database
        Cds: Cds database table
    """

    valid_cds_keys = {
        'id',
        'start',
        'end',
        'chromosome',
        'exon_id'
    }

    for cds, cds_items in gencode_cds_dict.items():
        filtered_cds_data = {}
        
        # Filtering the dict to only contain the key-value pairs that are contained in the database
        filtered_cds_data = {key: value for key, value in cds_items.items() if key in valid_cds_keys}

        # Using ** it recieves the following data: gene_id="gene_id", start=100, end=200 ...
        new_cds = Cds(**filtered_cds_data)

        try:
            print(f"inserting cds {new_cds}")
            # Add new gene to the session
            session.add(new_cds)
            session.commit()
        except:
            msg = "Error in database writing"
            raise SQLAlchemyError(msg)
        else:
            msg = f"INFO: gene {new_cds} record successfully"
            logging.info(msg)

# insert_db_gencode_gene_data(gencode_gene_dict, session, Genes)
# insert_refseq_transcripts_data(transcripts_ids, session, RefseqTranscripts)
# insert_db_gencode_transcirpts_data(gencode_transcript_dict, session, Transcripts, RefseqTranscripts, association_table_trans_coords, association_table_exons, association_table_cds)
# insert_db_gencode_exons_data(gencode_exon_dict, session, Exons)
# insert_db_gencode_cds_data(gencode_cds_dict, session, Cds)