from Bio import Entrez

def convert_refseq_to_ensembl(refseq_id):
    Entrez.email = "oriolcanal1998@gmail.com"  # Set your email address for NCBI API usage


    # Search the gene in NCBI Gene database using the RefSeq identifier
    handle = Entrez.esearch(db="gene", term=refseq_id)
    record = Entrez.read(handle)

    if record["IdList"]:
        gene_id = record["IdList"][0]

        # Retrieve the gene information from NCBI Gene database
        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        record = Entrez.read(handle)

        # Extract the Ensembl gene identifier if available
        for gene in record[0]["Entrezgene_gene"]:
            if isinstance(gene, dict) and gene.get("Gene-ref_db") == "Ensembl":
                ensembl_id = gene["Gene-ref_locus"]
                return ensembl_id

    return None  # Return None if conversion is not successful

# Example usage
refseq_gene_id = "NC_000001.11"
ensembl_gene_id = convert_refseq_to_ensembl(refseq_gene_id)

if ensembl_gene_id:
    print("Ensembl gene ID:", ensembl_gene_id)
else:
    print("No corresponding Ensembl gene ID found.")