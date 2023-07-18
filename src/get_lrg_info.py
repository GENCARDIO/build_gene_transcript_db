def get_lrg_trancripts(lrg_txt: str):
    lrg_ensembl = dict()
    lrg_refseq = dict()
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
            ccds = tmp[6].strip()

            if ensembl_trans_id not in lrg_ensembl:
                lrg_ensembl.setdefault(ensembl_trans_id, dict())
                lrg_ensembl[ensembl_trans_id]["lrg_id"] = lrg_id
                lrg_ensembl[ensembl_trans_id]["lrg_transcript"] = lrg_transcript
                lrg_ensembl[ensembl_trans_id]["ccds"] = ccds
            
            if refseq_transcript_id not in lrg_refseq:
                lrg_refseq.setdefault(refseq_transcript_id, dict())
                lrg_refseq[refseq_transcript_id]["lrg_id"] = lrg_id
                lrg_refseq[refseq_transcript_id]["lrg_transcript"] = lrg_transcript
                lrg_refseq[refseq_transcript_id]["ccds"] = ccds




    return (lrg_ensembl, lrg_refseq)

lrg_ensembl, lrg_refseq = get_lrg_trancripts("/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/lrg/list_LRGs_transcripts_xrefs.txt")


