from src.get_refseq_info import (
    refseq_genename_geneobject,
    refseq_genename_geneid_transobject,
    refseq_genename_geneid_transid_exonobject,
    refseq_genename_geneid_transid_cdsobject,
)
from src.get_gencode_info import (
    gencode_genename_geneobject,
    gencode_genename_geneid_transobject,
    gencode_genename_geneid_transid_exonobject,
    gencode_genename_geneid_transid_cdsobject
)
from src.parse_hugo import gene_name_to_synonym_dict
from src.global_variables import logging


def compare_refseq_ensembl_genes(
    gencode_genename_geneobject: dict,
    refseq_genename_geneobject: dict,
    gene_name_to_synonym_dict: dict
):
    logging.info("Comparing RefSeq and Ensembl genes")
    for gene_name, ensembl_gene_objects in gencode_genename_geneobject.items():
        gene_synonyms = set()
        gene_synonyms.add(gene_name)

        # geting the synonyms from the gene_name

        if gene_name in gene_name_to_synonym_dict:
            gene_synonyms.update(gene_name_to_synonym_dict[gene_name])

        # adding all gencode and refseq objects related to a gene_name or its synonyms
        ensembl_genes = list()
        refseq_genes = list()
        for gene_name in gene_synonyms:
            if gene_name in gencode_genename_geneobject:
                ensembl_genes.extend(gencode_genename_geneobject[gene_name])
            if gene_name in refseq_genename_geneobject:
                refseq_genes.extend(refseq_genename_geneobject[gene_name])
        
        # comparing coords of refseq and gencode 
        for ensembl_gene in ensembl_genes:
            ensembl_coords = ensembl_gene.get_coords()
            ensembl_gene.gene_synonyms = gene_synonyms
            for refseq_gene in refseq_genes:
                refseq_gene.gene_synonyms = gene_synonyms
                refseq_coords = refseq_gene.get_coords()
                if refseq_coords == ensembl_coords:
                    if ensembl_gene.refseq_same_gene_coords is None:
                        ensembl_gene.refseq_same_gene_coords = {refseq_gene.gene_id}
                    else:
                        ensembl_gene.refseq_same_gene_coords.add(refseq_gene.gene_id)
                    if refseq_gene.gencode_same_gene_coords is None:
                        refseq_gene.gencode_same_gene_coords = {ensembl_gene.gene_id}
                    else:
                        refseq_gene.gencode_same_gene_coords.add(ensembl_gene.gene_id)


def compare_refseq_ensembl_transcripts(
    gencode_genename_geneid_transobject: dict,
    refseq_genename_geneid_transobject: dict,
    gene_name_to_synonym_dict: dict
):
    logging.info("Comparing refseq and ensembl transcripts")
    for gene_name, ensembl_gene_objects in gencode_genename_geneid_transobject.items():
        gene_synonyms = set()
        gene_synonyms.add(gene_name)

        # geting the synonyms from the gene_name

        if gene_name in gene_name_to_synonym_dict:
            gene_synonyms.update(gene_name_to_synonym_dict[gene_name])

        # adding all gencode and refseq objects related to a gene_name or its synonyms
        ensembl_transcripts = list()
        refseq_transcripts = list()
        for gene_name in gene_synonyms:
            if gene_name in gencode_genename_geneid_transobject:
                for gencode_id in gencode_genename_geneid_transobject[gene_name]:
                    ensembl_transcripts.extend(gencode_genename_geneid_transobject[gene_name][gencode_id])
            if gene_name in refseq_genename_geneid_transobject:
                for refseq_id in refseq_genename_geneid_transobject[gene_name]:
                    refseq_transcripts.extend(refseq_genename_geneid_transobject[gene_name][refseq_id])

        # comparing coords of refseq and gencode 
        for ensembl_trans in ensembl_transcripts:
            ensembl_coords = ensembl_trans.get_coords()
            for refseq_trans in refseq_transcripts:
                refseq_coords = refseq_trans.get_coords()
                if refseq_coords == ensembl_coords:
                    if ensembl_trans.refseq_same_trans_coords is None:
                        ensembl_trans.refseq_same_trans_coords = {refseq_trans.id}
                    else:
                        ensembl_trans.refseq_same_trans_coords.add(refseq_trans.id)
                    if refseq_trans.gencode_same_trans_coords is None:
                        refseq_trans.gencode_same_trans_coords = {ensembl_trans.id}
                    else:
                        refseq_trans.gencode_same_trans_coords.add(ensembl_trans.id)


def compare_refseq_ensembl_exons(
    gencode_genename_geneid_transid_exonobject: dict,
    refseq_genename_geneid_transid_exonobject: dict,
    gene_name_to_synonym_dict: dict
):
    logging.info("Comparing RefSeq and Ensembl exons")

    for gene_name, _ in gencode_genename_geneid_transid_exonobject.items():
        gene_synonyms = set()
        gene_synonyms.add(gene_name)

        if gene_name in gene_name_to_synonym_dict:
            gene_synonyms.update(gene_name_to_synonym_dict[gene_name])

        # creating a dictionary to store transcript_id : [exon_objects]
        # associated to each transcript this dictionary will contain
        # transcripts and its associated exons that are referenced to
        #  a specific gene_name or its synonyms
        ensembl_transid_exons = dict()
        refseq_transid_exons = dict()
        for gene_name in gene_synonyms:
            if gene_name in gencode_genename_geneid_transid_exonobject:
                for gencode_id in gencode_genename_geneid_transid_exonobject[gene_name]:
                    for transcript_obj in gencode_genename_geneid_transid_exonobject[gene_name][gencode_id]:
                        for exon_objec in gencode_genename_geneid_transid_exonobject[gene_name][gencode_id][transcript_obj]:
                            if transcript_obj in ensembl_transid_exons:
                                ensembl_transid_exons[transcript_obj].append(exon_objec)
                            else:
                                ensembl_transid_exons[transcript_obj] = [exon_objec]
            if gene_name in refseq_genename_geneid_transid_exonobject:
                for refseq_id in refseq_genename_geneid_transid_exonobject[gene_name]:
                    for refseq_trans_obj in refseq_genename_geneid_transid_exonobject[gene_name][refseq_id]:
                        for refseq_exon_obj in refseq_genename_geneid_transid_exonobject[gene_name][refseq_id][refseq_trans_obj]:
                            if refseq_trans_obj in refseq_transid_exons:
                                refseq_transid_exons[refseq_trans_obj].append(refseq_exon_obj)
                            else:
                                refseq_transid_exons[refseq_trans_obj] = [refseq_exon_obj]

        # comparing coords of refseq and gencode 
        for ensembl_trans_obj, ensembl_exons in ensembl_transid_exons.items():
            ensembl_transcript_exon_coords = set()
            # get all exon coords from a transcript
            for ensembl_exon_obj in ensembl_exons:
                ensembl_transcript_exon_coords.add(ensembl_exon_obj.get_coords()) 

                for refseq_trans_obj, refseq_exons in refseq_transid_exons.items():
                    refseq_transcript_exon_coords = set()
                    for refseq_exon_obj in refseq_exons:
                        refseq_transcript_exon_coords.add(refseq_exon_obj.get_coords())
                    if refseq_transcript_exon_coords == ensembl_transcript_exon_coords:
                        # if the transcript class instance has refseq-same_exons_coords = None (default value),
                        # the value will be change to a list containing the refseq_gene id that matches the coordinates
                        if ensembl_trans_obj.refseq_same_exons_coords is not None:
                            ensembl_trans_obj.refseq_same_exons_coords.append(refseq_trans_obj.id)
                        else:
                            ensembl_trans_obj.refseq_same_exons_coords = [refseq_trans_obj.id]


def compare_refseq_ensembl_cds(
    gencode_genename_geneid_transid_cdsobject: dict,
    refseq_genename_geneid_transid_cdsobject: dict,
    gene_name_to_synonym_dict: dict
):
    logging.info("Comparing RefSeq and Ensembl cds")

    for gene_name, ensembl_gene in gencode_genename_geneid_transid_cdsobject.items():
        gene_synonyms = set()
        gene_synonyms.add(gene_name)

        # geting the synonyms from the gene_name

        if gene_name in gene_name_to_synonym_dict:
            gene_synonyms.update(gene_name_to_synonym_dict[gene_name])

        # creating a dictionary to store transcript_id : [exon_objects]
        # associated to each transcript this dictionary will contain
        # transcripts and its associated exons that are referenced to a
        # specific gene_name or its synonyms
        ensembl_transid_cds = dict()
        refseq_transid_cds = dict()
        for gene_name in gene_synonyms:
            if gene_name in gencode_genename_geneid_transid_cdsobject:
                for gencode_id in gencode_genename_geneid_transid_cdsobject[gene_name]:
                    for transcript_obj in gencode_genename_geneid_transid_cdsobject[gene_name][gencode_id]:
                        for cds_objec in gencode_genename_geneid_transid_cdsobject[gene_name][gencode_id][transcript_obj]:
                            if transcript_obj in ensembl_transid_cds:
                                ensembl_transid_cds[transcript_obj].append(cds_objec)
                            else:
                                ensembl_transid_cds[transcript_obj] = [cds_objec]
            if gene_name in refseq_genename_geneid_transid_cdsobject:
                for refseq_id in refseq_genename_geneid_transid_cdsobject[gene_name]:
                    for refseq_trans_obj in refseq_genename_geneid_transid_cdsobject[gene_name][refseq_id]:
                        for refseq_cds_obj in refseq_genename_geneid_transid_cdsobject[gene_name][refseq_id][refseq_trans_obj]:
                            if refseq_trans_obj in refseq_transid_cds:
                                refseq_transid_cds[refseq_trans_obj].append(refseq_cds_obj)
                            else:
                                refseq_transid_cds[refseq_trans_obj] = [refseq_cds_obj]

        # comparing coords of refseq and gencode 
        for ensembl_trans_obj, ensembl_cds in ensembl_transid_cds.items():
            ensembl_transcript_cds_coords = set()
            # get all exon coords from a transcript
            for ensembl_cds_obj in ensembl_cds:
                ensembl_transcript_cds_coords.add(ensembl_cds_obj.get_coords()) 
                
                for refseq_trans_obj, refseq_cds in refseq_transid_cds.items():
                    refseq_transcript_cds_coords = set()
                    for refseq_cds_obj in refseq_cds:
                        refseq_transcript_cds_coords.add(refseq_cds_obj.get_coords())
                    if refseq_transcript_cds_coords == ensembl_transcript_cds_coords:
                        # if the transcript class instance has refseq-same_exons_coords = None (default value),
                        # the value will be change to a list containing the refseq_gene id that matches the coordinates
                        if ensembl_trans_obj.refseq_same_cds_coords is not None:
                            ensembl_trans_obj.refseq_same_cds_coords.append(refseq_trans_obj.id)

                        else:
                            ensembl_trans_obj.refseq_same_cds_coords = [refseq_trans_obj.id]

                        
            # print(ensembl_trans_obj.refseq_same_exons_coords, "same_exons_coords")


# for gene_name, dict1 in refseq_genename_geneid_transid_exonobject.items():
#     for gene_id, dict2 in dict1.items():
#         for trans_obj, exons in dict2.items():
#             if trans_obj.transcript_id == "NM_153437.3":
#                 print("found")
#                 print(trans_obj.get_coords())

compare_refseq_ensembl_cds(
    gencode_genename_geneid_transid_cdsobject,
    refseq_genename_geneid_transid_cdsobject,
    gene_name_to_synonym_dict
)

compare_refseq_ensembl_exons(
    gencode_genename_geneid_transid_exonobject,
    refseq_genename_geneid_transid_exonobject,
    gene_name_to_synonym_dict
)

compare_refseq_ensembl_transcripts(
    gencode_genename_geneid_transobject,
    refseq_genename_geneid_transobject,
    gene_name_to_synonym_dict
)

compare_refseq_ensembl_genes(
    gencode_genename_geneobject,
    refseq_genename_geneobject,
    gene_name_to_synonym_dict
)
