import csv
from src.get_refseq_info import gene_name_to_synonym_dict
from global_variables import logging
hugo_path = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/hugo/hgnc_complete_set.txt"


def obtain_gene_synonyms_from_hugo(hugo_tsv, gene_name_to_synonym_dict):
    logging.info(f"importing gene synonyms from HUGO file: {hugo_tsv}")
    with open(hugo_tsv, "r") as file:
        reader = csv.reader(file, delimiter="\t")

        field_names = next(reader)

        for row in reader:
            all_synonyms = set()
            line_dict = dict(zip(field_names, row))
            all_synonyms.add(line_dict["symbol"])
            if "|" in line_dict["alias_symbol"]:
                alias_symbols = line_dict["alias_symbol"].split("|")
                all_synonyms.update(alias_symbols)
            else:
                all_synonyms.add(line_dict["alias_symbol"])
            if "|" in line_dict["prev_symbol"]:
                prev_symbol = line_dict["prev_symbol"].split("|")
                all_synonyms.update(prev_symbol)
            else:
                all_synonyms.add(line_dict["prev_symbol"])

            # removing empty elements from the set. e.g. if it has no
            # prev_symbol, "" will be added to the set and we want to remove it
            all_synonyms = {element for element in all_synonyms if element}
            for gene_sym in all_synonyms:
                if gene_sym not in gene_name_to_synonym_dict:
                    gene_name_to_synonym_dict[gene_sym] = all_synonyms
                else:
                    gene_name_to_synonym_dict[gene_sym].update(all_synonyms)
    logging.info(f"Gene synonyms have been obtained successfully!")
    return (gene_name_to_synonym_dict)


gene_name_to_synonym_dict = obtain_gene_synonyms_from_hugo(
    hugo_path,
    gene_name_to_synonym_dict
)
