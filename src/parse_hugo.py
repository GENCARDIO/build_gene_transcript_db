import csv
from src.get_refseq_info import parse_refseq_gff, refseq_gff_grch38
from src.get_mane import get_ensembl_mane, mane_ensembl_gff
from src.get_lrg_info import get_lrg_trancripts, lrg_gff
from src.global_variables import logging

hugo_path = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/hugo/hgnc_complete_set.txt"


def get_hugo_conf(
    hugo_path: str,
    versions_dict: dict
):
    versions_dict["resources"]["hugo"]["path"] = hugo_path
    return versions_dict


def obtain_gene_synonyms_from_hugo(hugo_tsv, gene_name_to_synonym_dict):
    logging.info(f"importing gene synonyms from HUGO file: {hugo_tsv}")
    with open(hugo_tsv, "r") as file:
        reader = csv.reader(file, delimiter="\t")

        field_names = next(reader)

        for row in reader:
            all_synonyms = set()
            line_dict = dict(zip(field_names, row))
            print(line_dict["symbol"], "symbol")

            all_synonyms.add(line_dict["symbol"])
            if "|" in line_dict["alias_symbol"]:
                alias_symbols = line_dict["alias_symbol"].split("|")
                all_synonyms.update(alias_symbols)
            else:
                all_synonyms.add(line_dict["alias_symbol"].upper())
            if "|" in line_dict["prev_symbol"]:
                prev_symbol = line_dict["prev_symbol"].split("|")
                all_synonyms.update(prev_symbol)
            else:
                all_synonyms.add(line_dict["prev_symbol"])

            # removing empty elements from the set. e.g. if it has no
            # prev_symbol, "" will be added to the set and we want to remove it
            all_synonyms = {element.upper() for element in all_synonyms if element}
            print(all_synonyms)
            for gene_sym in all_synonyms:
                if gene_sym not in gene_name_to_synonym_dict:
                    gene_name_to_synonym_dict[gene_sym] = all_synonyms
                else:
                    gene_name_to_synonym_dict[gene_sym].update(all_synonyms)
    logging.info(
        f"Gene synonyms from {hugo_path}\
        have been obtained successfully!")
    return (gene_name_to_synonym_dict)

if "__main__" == __name__:
    lrg_ensembl, lrg_refseq = get_lrg_trancripts(lrg_gff)

    (
        mane_clin_ens_id,
        mane_select_ens_id,
        mane_clin_refseq_id,
        mane_select_refseq_id
    ) = get_mane(mane_gff)

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
    gene_name_to_synonym_dict = obtain_gene_synonyms_from_hugo(
        hugo_path,
        gene_name_to_synonym_dict
    )