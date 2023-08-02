from src.get_refseq_info import (
    parse_refseq_gff,
    refseq_gff_grch38,
    refseq_gff_grch37,
    mock_refseq_gff_grch38,
    mock_refseq_gff_grch37
)
from src.get_mane import get_mane, mane_gff
from src.get_gencode_info import (
    parse_gencode,
    gencode_gff_grch38,
    gencode_gff_grch37,
    mock_gencode_gff_grch37,
    mock_gencode_gff_grch38
)
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
from src.parse_hugo import (
    obtain_gene_synonyms_from_hugo,
    hugo_path
)
from src.create_database import create_mysql_database

if "__main__" == __name__:
    gene_name_to_synonym_dict = dict()
    genome_v = 38
    # extracting refseq and ensembl ids that are mane clin or select
    (
        mane_clin_ens_id,
        mane_select_ens_id,
        mane_clin_refseq_id,
        mane_select_refseq_id
    ) = get_mane(mane_gff)

    # extracting lrg_id, lrg_transcript and ccds
    lrg_ensembl, lrg_refseq = get_lrg_trancripts(lrg_gff)
    
    # extracting refseq elements
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
        lrg_refseq,
        genome_v,
        gene_name_to_synonym_dict
    )

    #
    gene_name_to_synonym_dict = obtain_gene_synonyms_from_hugo(
        hugo_path,
        gene_name_to_synonym_dict
    )

    # extracting gencode elements
    (
        gencode_genename_geneobject,
        gencode_genename_geneid_transobject,
        gencode_genename_geneid_transid_exonobject,
        gencode_genename_geneid_transid_cdsobject
    ) = parse_gencode(
        gencode_gff_grch38,
        mane_clin_ens_id,
        mane_select_ens_id,
        lrg_ensembl,
        genome_v
    )
    
    # GENOME VERSION = 37 -------------------------------------
    genome_v = 37

    # extracting refseq elements
    (
        refseq_genename_geneobject_grch37,
        refseq_genename_geneid_transobject_grch37,
        refseq_genename_geneid_transid_exonobject_grch37,
        refseq_genename_geneid_transid_cdsobject_grch37,
        gene_name_to_synonym_dict
    ) = parse_refseq_gff(
        refseq_gff_grch37,
        mane_clin_refseq_id,
        mane_select_refseq_id,
        lrg_refseq,
        genome_v,
        gene_name_to_synonym_dict
    )

    # extracting gencode elements
    (
        gencode_genename_geneobject_grch37,
        gencode_genename_geneid_transobject_grch37,
        gencode_genename_geneid_transid_exonobject_grch37,
        gencode_genename_geneid_transid_cdsobject_grch37
    ) = parse_gencode(
        gencode_gff_grch37,
        mane_clin_ens_id,
        mane_select_ens_id,
        lrg_ensembl,
        genome_v
    )

    # --------------------------COMPARING REFSEQ-GENCODE----------------------

    # comparing refseq and gencode elements GRCH38:
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


    # comparing refseq and gencode elements GRCH37:
    compare_refseq_ensembl_cds(
        gencode_genename_geneid_transid_cdsobject_grch37,
        refseq_genename_geneid_transid_cdsobject_grch37,
        gene_name_to_synonym_dict
    )

    compare_refseq_ensembl_exons(
        gencode_genename_geneid_transid_exonobject_grch37,
        refseq_genename_geneid_transid_exonobject_grch37,
        gene_name_to_synonym_dict
    )

    compare_refseq_ensembl_transcripts(
        gencode_genename_geneid_transobject_grch37,
        refseq_genename_geneid_transobject_grch37,
        gene_name_to_synonym_dict
    )

    compare_refseq_ensembl_genes(
        gencode_genename_geneobject_grch37,
        refseq_genename_geneobject_grch37,
        gene_name_to_synonym_dict
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
        # RefseqGenes_Grch37,
        # RefseqTranscripts_Grch37,
        # RefseqExons_Grch37,
        # RefseqCds_Grch37,
        # Genes_Grch37,
        # Transcripts_Grch37,
        # Exons_Grch37,
        # Cds_Grch37,
        gencode_gene_synonyms_association,
        association_table_trans_coords,
        association_table_exon_coords,
        association_table_cds_coords,
        session
    ) = create_mysql_database()

    insert_db_refseq_genes(  # grch38
        refseq_genename_geneobject,
        session,
        RefseqGenes,
        Gene_Synonyms
    )
    insert_db_refseq_genes(  # grch37
        refseq_genename_geneobject_grch37,
        session,
        RefseqGenes,
        Gene_Synonyms
    )

    insert_db_gencode_gene_data(  # grch38
        gencode_genename_geneobject,
        session,
        Genes,
        Gene_Synonyms,
        gencode_gene_synonyms_association,
        RefseqGenes
    )
    insert_db_gencode_gene_data(  # grch37
        gencode_genename_geneobject_grch37,
        session,
        Genes,
        Gene_Synonyms,
        gencode_gene_synonyms_association,
        RefseqGenes
    )

    insert_db_refseq_transcripts(  # grch38
        refseq_genename_geneid_transobject,
        session,
        RefseqTranscripts
    )
    insert_db_refseq_transcripts(  # grch37
        refseq_genename_geneid_transobject_grch37,
        session,
        RefseqTranscripts
    )
    insert_db_refseq_exons(  # grch38
        refseq_genename_geneid_transid_exonobject,
        session,
        RefseqExons
    )
    insert_db_refseq_exons(  # grch37
        refseq_genename_geneid_transid_exonobject_grch37,
        session,
        RefseqExons
    )
    insert_db_refseq_cds(  # grch38
        refseq_genename_geneid_transid_cdsobject,
        session,
        RefseqCds
    )
    insert_db_refseq_cds(  # grch37
        refseq_genename_geneid_transid_cdsobject_grch37,
        session,
        RefseqCds
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
        RefseqTranscripts
    )
    insert_db_gencode_exons(  #grch38
        gencode_genename_geneid_transid_exonobject,
        session,
        Exons
    )
    insert_db_gencode_exons(  #grch37
        gencode_genename_geneid_transid_exonobject_grch37,
        session,
        Exons
    )
    insert_db_gencode_cds(  # grch38
        gencode_genename_geneid_transid_cdsobject,
        session,
        Cds
    )
    insert_db_gencode_cds(  # grch37
        gencode_genename_geneid_transid_cdsobject_grch37,
        session,
        Cds
    )
