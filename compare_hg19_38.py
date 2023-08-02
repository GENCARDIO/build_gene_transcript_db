from src.get_refseq_info import (
    parse_refseq_gff,
    refseq_gff_grch38,
    refseq_gff_grch37
)
from src.get_mane import get_mane, mane_gff
from src.get_gencode_info import parse_gencode, gencode_gff
from src.get_lrg_info import get_lrg_trancripts, lrg_gff
from src.compare_gencode_refseq import (
    compare_refseq_ensembl_cds,
    compare_refseq_ensembl_exons,
    compare_refseq_ensembl_transcripts,
    compare_refseq_ensembl_genes
)
from src.inserting_into_db import (
    insert_db_gencode_gene_data,
    insert_db_refseq_genes,
    insert_db_refseq_transcripts,
    insert_db_refseq_exons,
    insert_db_refseq_cds,
    insert_db_gencode_transcripts,
    insert_db_gencode_exons,
    insert_db_gencode_cds
)
from src.create_database import create_mysql_database

if "__main__" == __name__:
    # extracting refseq and ensembl ids that are mane clin or select
    (
        mane_clin_ens_id,
        mane_select_ens_id,
        mane_clin_refseq_id,
        mane_select_refseq_id
    ) = get_mane(mane_gff)

    # extracting lrg_id, lrg_transcript and ccds
    lrg_ensembl, lrg_refseq = get_lrg_trancripts(lrg_gff)

    # extracting refseq grch38 elements
    (
        refseq_genename_geneobject,
        refseq_genename_geneid_transobject,
        refseq_genename_geneid_transid_exonobject,
        refseq_genename_geneid_transid_cdsobject,
        gene_name_to_synonym_dict
    ) = parse_refseq_gff(
        refseq_gff_grch38,
        mane_clin_refseq_id,
        mane_select_refseq_id,
        lrg_refseq
    )

    # extracting refseq grch37 elements
    (
        refseq_genename_geneobject_grch37,
        refseq_genename_geneid_transobject_grch37,
        refseq_genename_geneid_transid_exonobject_grch37,
        refseq_genename_geneid_transid_cdsobject_grch37,
        gene_name_to_synonym_dict_grch37
    ) = parse_refseq_gff(
        refseq_gff_grch37,
        mane_clin_refseq_id,
        mane_select_refseq_id,
        lrg_refseq
    )
    # print(gene_name_to_synonym_dict, "hg38")
    # print(gene_name_to_synonym_dict_grch37["MIR1302-2HG"], "hg37")
    grch37_objs = list()
    # for gene_name, gene_objs in refseq_genename_geneobject_grch37.items():
    #     for gene_obj in gene_objs:
    #         grch37_objs.append(gene_obj)

    # for gene_name, gene_objs2 in refseq_genename_geneobject.items():
    #     for gene_obj2 in gene_objs2:
    #         gene_id = gene_obj2.gene_id
    #         gene_id_found = False
    #         for gene_obj in grch37_objs:
    #             if gene_obj.gene_id == gene_id:
    #                 gene_id_found = True
    #         print(gene_id_found, gene_id)


    for gene_name, dict in refseq_genename_geneid_transobject_grch37.items():
        for gene_id, transobjs in dict.items():
            for trans_obj in transobjs:
                grch37_objs.append(trans_obj)

    for gene_name, dict2 in refseq_genename_geneid_transobject.items():
        for gene_id, transobjs2 in dict2.items():
            for trans_obj2 in transobjs2:
                found = False
                trans_id = trans_obj2.id.rsplit(".", 1)[0]

                for trans_obj in grch37_objs:
                    if trans_obj.id.rsplit(".", 1)[0] == trans_id:
                        found = True
                print(found, trans_id)