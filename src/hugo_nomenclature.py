# import requests

# hgnc_id = "HGNC:1100"

# server = "https://rest.ensembl.org"
# ext = f"/xrefs/id/{hgnc_id}?object_type=gene"
# headers = {"Content-Type": "application/json"}
# response = requests.get(server + ext, headers=headers)
# print(response)
# if response.ok:
#     data = response.json()
#     ensembl_gene_id = None
#     for entry in data:
#         if entry["dbname"] == "Ensembl Gene":
#             ensembl_gene_id = entry["id"]
#             break

#     if ensembl_gene_id:
#         print(f"Ensembl Gene ID for {hgnc_id}: {ensembl_gene_id}")
#     else:
#         print(f"No Ensembl Gene ID found for {hgnc_id}")
# else:
#     print(f"Request failed with status code {response.status_code}")

# from PyEntrezId import Conversion

# EnsemblId = "ENST00000407559"
# Id = Conversion("ocanal@idibgi.org")

# EntrezId = Id.convert_ensembl_to_entrez(EnsemblId)

# print(EntrezId)

import sys
import mygene

# mg = mygene.MyGeneInfo()

# genes = ["MIR6859-1", "MIR1302-2HG"]


# for gene in genes:
#     result = mg.query(gene, scopes="symbol", fields=["ensembl"], species="human", verbose=False)
#     hgnc_name = gene
#     for hit in result["hits"]:
#         if "ensembl" in hit and "gene" in hit["ensembl"]:
#             sys.stdout.write("%s\t%s\n" % (hgnc_name, hit["ensembl"]["gene"]))
mg = mygene.MyGeneInfo()

g = mg.getgene(102466751)# query by Entrez id
print(g[ 'entrezgene']) # same as input 
print()
print( g["genomic_pos"]["ensemblgene"] )