import datetime
import yaml
import os

from src.global_variables import logging, get_conf
from src.get_refseq_info import (
    parse_refseq_gff,
    get_refseq_conf,
)
from src.get_mane import get_ensembl_mane, get_refseq_mane, get_mane_conf
from src.get_gencode_info import (
    parse_gencode,
    get_gencode_conf,
)
from src.get_lrg_info import get_lrg_trancripts, get_lrg_conf
from src.compare_gencode_refseq import (
    compare_refseq_ensembl_cds,
    compare_refseq_ensembl_exons,
    compare_refseq_ensembl_transcripts,
    compare_refseq_ensembl_genes,
)
from src.inserting_into_db import (
    insert_db_gencode_gene_data,
    insert_db_refseq_genes,
    insert_db_refseq_transcripts,
    insert_db_refseq_exons,
    insert_db_refseq_cds,
    insert_db_gencode_transcripts,
    insert_db_gencode_exons,
    insert_db_gencode_cds,
)
from src.parse_hugo import obtain_gene_synonyms_from_hugo, get_hugo_conf
from src.create_database import create_mysql_database


def create_conf(
    refseq_gff_grch38,
    refseq_gff_grch37,
    gencode_gff_grch38,
    gencode_gff_grch37,
    mane_refseq_gff,
    mane_ensembl_gff,
    hugo_path,
    lrg_gff,
    main_path,
):
    datetime_obj = datetime.datetime.now()
    current_datetime = datetime_obj.strftime("%Y-%m-%d:%Hh:%Mm:%Ss")
    versions_dict = {
        "database_name": genes_db_name,
        "date_created": current_datetime,
        "resources": {
            "refseq": {
                "grch37": {
                    "genome_build": None,
                    "accession_version": None,
                    "file_path": None,
                },
                "grch38": {
                    "genome_build": None,
                    "accession_version": None,
                    "file_path": None,
                },
            },
            "gencode": {
                "grch37": {
                    "gencode_version": None,
                    "ensembl_version": None,
                    "file_path": None,
                },
                "grch38": {
                    "gencode_version": None,
                    "ensembl_version": None,
                    "file_path": None,
                },
            },
            "mane": {
                "grch37": {
                    "refseq": {
                        "version": None,
                        "file_path": None,
                    },
                    "ensembl": {"version": None, "file_path": None},
                },
                "grch38": {
                    "refseq": {
                        "version": None,
                        "file_path": None,
                    },
                    "ensembl": {"version": None, "file_path": None},
                },
            },
            "hugo": {"file_path": None},
            "lrg": {"file_path": None, "version": None},
        },
    }
    versions_dict = get_refseq_conf(refseq_gff_grch38, versions_dict)
    versions_dict = get_refseq_conf(refseq_gff_grch37, versions_dict)
    versions_dict = get_gencode_conf(gencode_gff_grch38, versions_dict)
    versions_dict = get_gencode_conf(gencode_gff_grch37, versions_dict)
    versions_dict = get_mane_conf(mane_refseq_gff, versions_dict, 38)
    versions_dict = get_mane_conf(mane_ensembl_gff, versions_dict, 38)
    versions_dict = get_hugo_conf(hugo_path, versions_dict)
    versions_dict = get_lrg_conf(lrg_gff, versions_dict)

    conf_filename = f"{genes_db_name}_database_conf.yaml"
    conf_file_path = os.path.join(main_path, conf_filename)
    print(conf_file_path)
    with open(conf_file_path, "w") as file:
        yaml.dump(versions_dict, file)

    logging.info(f"config file created in {conf_file_path}")


if "__main__" == __name__:
    # conf file where paths are stored
    conf_file = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/resources_config.yaml"

    (
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
        main_path,
    ) = get_conf(conf_file)

    create_conf(
        refseq_gff_grch38,
        refseq_gff_grch37,
        gencode_gff_grch38,
        gencode_gff_grch37,
        mane_refseq_gff,
        mane_ensembl_gff,
        hugo_path,
        lrg_gff,
        main_path,
    )

    gene_name_to_synonym_dict = dict()
    files_version = dict()
    genome_v = 38
    # extracting refseq and ensembl ids that are mane clin or select
    (
        mane_clin_ensembl_id,
        mane_clin_refseq_id,
        mane_select_ensembl_refseq,
        mane_select_refseq_ensembl,
    ) = get_ensembl_mane(mane_ensembl_gff)
    # print(len(mane_select_ensembl_refseq))

    (
        mane_clin_ensembl_id,
        mane_clin_refseq_id,
        mane_select_ensembl_refseq,
        mane_select_refseq_ensembl,
    ) = get_refseq_mane(
        mane_refseq_gff,
        mane_clin_refseq_id,
        mane_clin_ensembl_id,
        mane_select_refseq_ensembl,
        mane_select_ensembl_refseq,
    )

    # extracting lrg_id, lrg_transcript and ccds
    lrg_ensembl, lrg_refseq = get_lrg_trancripts(lrg_gff)

    # extracting refseq elements
    (
        refseq_genename_geneobject,
        refseq_genename_geneid_transobject,
        refseq_genename_geneid_transid_exonobject,
        refseq_genename_geneid_transid_cdsobject,
        gene_name_to_synonym_dict,
    ) = parse_refseq_gff(
        refseq_gff_grch38,
        mane_clin_refseq_id,
        mane_select_refseq_ensembl,
        lrg_refseq,
        genome_v,
        gene_name_to_synonym_dict,
    )

    #

    # extracting gencode elements
    (
        gencode_genename_geneobject,
        gencode_genename_geneid_transobject,
        gencode_genename_geneid_transid_exonobject,
        gencode_genename_geneid_transid_cdsobject,
    ) = parse_gencode(
        gencode_gff_grch38,
        mane_clin_ensembl_id,
        mane_select_ensembl_refseq,
        lrg_ensembl,
        genome_v,
    )

    # # GENOME VERSION = 37 -------------------------------------
    genome_v = 37

    # extracting refseq elements
    (
        refseq_genename_geneobject_grch37,
        refseq_genename_geneid_transobject_grch37,
        refseq_genename_geneid_transid_exonobject_grch37,
        refseq_genename_geneid_transid_cdsobject_grch37,
        gene_name_to_synonym_dict,
    ) = parse_refseq_gff(
        refseq_gff_grch37,
        mane_clin_refseq_id,
        mane_select_refseq_ensembl,
        lrg_refseq,
        genome_v,
        gene_name_to_synonym_dict,
    )

    # # extracting gencode elements
    (
        gencode_genename_geneobject_grch37,
        gencode_genename_geneid_transobject_grch37,
        gencode_genename_geneid_transid_exonobject_grch37,
        gencode_genename_geneid_transid_cdsobject_grch37,
    ) = parse_gencode(
        gencode_gff_grch37,
        mane_clin_ensembl_id,
        mane_select_ensembl_refseq,
        lrg_ensembl,
        genome_v,
    )

    gene_name_to_synonym_dict = obtain_gene_synonyms_from_hugo(
        hugo_path, gene_name_to_synonym_dict
    )

    print(gene_name_to_synonym_dict["TRNL2"])

    # --------------------------COMPARING REFSEQ-GENCODE----------------------

    # comparing refseq and gencode elements GRCH38:
    compare_refseq_ensembl_cds(
        gencode_genename_geneid_transid_cdsobject,
        refseq_genename_geneid_transid_cdsobject,
        gene_name_to_synonym_dict,
    )

    compare_refseq_ensembl_exons(
        gencode_genename_geneid_transid_exonobject,
        refseq_genename_geneid_transid_exonobject,
        gene_name_to_synonym_dict,
    )

    compare_refseq_ensembl_transcripts(
        gencode_genename_geneid_transobject,
        refseq_genename_geneid_transobject,
        gene_name_to_synonym_dict,
    )

    compare_refseq_ensembl_genes(
        gencode_genename_geneobject,
        refseq_genename_geneobject,
        gene_name_to_synonym_dict,
    )

    # comparing refseq and gencode elements GRCH37:
    compare_refseq_ensembl_cds(
        gencode_genename_geneid_transid_cdsobject_grch37,
        refseq_genename_geneid_transid_cdsobject_grch37,
        gene_name_to_synonym_dict,
    )

    compare_refseq_ensembl_exons(
        gencode_genename_geneid_transid_exonobject_grch37,
        refseq_genename_geneid_transid_exonobject_grch37,
        gene_name_to_synonym_dict,
    )

    compare_refseq_ensembl_transcripts(
        gencode_genename_geneid_transobject_grch37,
        refseq_genename_geneid_transobject_grch37,
        gene_name_to_synonym_dict,
    )

    compare_refseq_ensembl_genes(
        gencode_genename_geneobject_grch37,
        refseq_genename_geneobject_grch37,
        gene_name_to_synonym_dict,
    )

    # create database:
    (
        Gene_Synonyms,
        RefseqGenes,
        RefseqTranscripts,
        RefseqExons,
        RefseqCds,
        Genes,
        Transcripts,
        Exons,
        Cds,
        File_version,
        # RefseqGenes_Grch37,
        # RefseqTranscripts_Grch37,
        # RefseqExons_Grch37,
        # RefseqCds_Grch37,
        # Genes_Grch37,
        # Transcripts_Grch37,
        # Exons_Grch37,
        # Cds_Grch37,
        refseq_gene_synonyms_association,
        gencode_gene_synonyms_association,
        association_table_trans_coords,
        association_table_exon_coords,
        association_table_cds_coords,
        session,
    ) = create_mysql_database()

    insert_db_refseq_genes(  # grch38
        refseq_genename_geneobject, session, RefseqGenes, Gene_Synonyms
    )
    insert_db_refseq_genes(  # grch37
        refseq_genename_geneobject_grch37, session, RefseqGenes, Gene_Synonyms
    )

    insert_db_gencode_gene_data(  # grch38
        gencode_genename_geneobject,
        session,
        Genes,
        Gene_Synonyms,
        gencode_gene_synonyms_association,
        RefseqGenes,
    )
    insert_db_gencode_gene_data(  # grch37
        gencode_genename_geneobject_grch37,
        session,
        Genes,
        Gene_Synonyms,
        gencode_gene_synonyms_association,
        RefseqGenes,
    )

    insert_db_refseq_transcripts(  # grch38
        refseq_genename_geneid_transobject, session, RefseqTranscripts
    )
    insert_db_refseq_transcripts(  # grch37
        refseq_genename_geneid_transobject_grch37, session, RefseqTranscripts
    )
    insert_db_refseq_exons(  # grch38
        refseq_genename_geneid_transid_exonobject, session, RefseqExons
    )
    insert_db_refseq_exons(  # grch37
        refseq_genename_geneid_transid_exonobject_grch37, session, RefseqExons
    )
    insert_db_refseq_cds(  # grch38
        refseq_genename_geneid_transid_cdsobject, session, RefseqCds
    )
    insert_db_refseq_cds(  # grch37
        refseq_genename_geneid_transid_cdsobject_grch37, session, RefseqCds
    )
    insert_db_gencode_transcripts(  # grch38
        gencode_genename_geneid_transobject,
        session,
        Transcripts,
        RefseqTranscripts
    )
    insert_db_gencode_transcripts(  # grch37
        gencode_genename_geneid_transobject_grch37,
        session,
        Transcripts,
        RefseqTranscripts,
    )
    insert_db_gencode_exons(  # grch38
        gencode_genename_geneid_transid_exonobject, session, Exons
    )
    insert_db_gencode_exons(  # grch37
        gencode_genename_geneid_transid_exonobject_grch37, session, Exons
    )
    insert_db_gencode_cds(  # grch38
        gencode_genename_geneid_transid_cdsobject, session, Cds
    )
    insert_db_gencode_cds(  # grch37
        gencode_genename_geneid_transid_cdsobject_grch37, session, Cds
    )
