from dataclasses import dataclass, field
import gzip


from src.global_variables import convert_chromosomes, logging
from src.get_mane import mane_clin_refseq_id, mane_select_refseq_id
from src.get_lrg_info import lrg_refseq

refseq_gff = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/refseq/GRCh38_latest_genomic.gff.gz"


@dataclass
class RefseqGene:
    gene_id:                    str
    start:                      int
    end:                        int
    chromosome:                 int
    id:                         str
    name:                       str
    gene_biotype:               str
    feature:                    str
    gene_synonyms:              str = field(default=None)
    hgnc_id:                    int = field(default=None)
    num_transcripts:            int = field(default=0)
    gencode_same_gene_coords:   str = field(default=None)

    def get_coords(self):
        return f"{self.chromosome}:{self.start}-{self.end}"

    def increment_num_transcripts(self):
        self.num_transcripts += 1


@dataclass
class RefseqTranscript:
    id:                         str
    start:                      int
    end:                        int
    chromosome:                 str
    gene_id:                    int
    gene_name:                  str
    hgnc_id:                    str = field(default=None)
    lrg_id:                     str = field(default=None)
    lrg_transcript:             str = field(default=None)
    ccds:                       str = field(default=None)
    numb_exons:                 int = field(default=0)
    numb_cds:                   int = field(default=0)
    mane_clin:                  bool = field(default=False)
    mane_select:                bool = field(default=False)
    gencode_same_trans_coords:  str = field(default=None)
    gencode_same_cds_coords:    str = field(default=None)

    def increment_numb_exons(self):
        self.numb_exons += 1

    def increment_numb_cds(self):
        self.numb_cds += 1

    def get_coords(self):
        return f"{self.chromosome}:{self.start}-{self.end}"

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        if isinstance(other, RefseqTranscript):
            return self.id == other.id
        return False


@dataclass
class RefseqExons:
    exon_id:        str
    transcript_id:  str
    start:          int
    end:            int
    chromosome:     str
    exon_number:    int

    def get_coords(self):
        return f"{self.chromosome}:{self.start}-{self.end}"


@dataclass
class RefseqCds:
    id:             str
    start:          int
    end:            int
    chromosome:     str
    transcript_id:  str

    def get_coords(self):
        return f"{self.chromosome}:{self.start}-{self.end}"

def parse_gff3_common_fields(line):
    feature_dict = dict()
    line = line.rstrip("\n")
    if line.startswith("#"):
        continue
    tmp = line.split("\t")
    feature_dict["feature"] = tmp[2].lower()
    refseq_chr = tmp[0]
    feature_dict["chromosome"] = convert_chromosomes(refseq_chr)
    feature_dict["start"] = tmp[3]
    feature_dict["end"] = tmp[4]
    feature_dict["info"] = tmp[8].split(";")
    feature_dict["chromosome"] = chr
    feature_dict["start"] = pos
    feature_dict["end"] = end

    return (feature_dict)
def get_refseq_transcripts(
    refseq_gff: str,
    mane_clin_refseq_id: list,
    mane_select_refseq_id: list
):
    '''
        Function that reads a RefSeq gff3 file and returns 5 dictionaries
        to identify genes, transcripts, exons and cds associated with the
        same gene name:

        params:
        -------
            refseq_gff: path of the refseq gff file
            mane_clin_refseq_id: list containing transcript refseq ids that
                are identified as mane clin in the mane gff3 file
            mane_select_refseq_id: list containing transcript refseq ids that
                are identified as mane select in the mane gff3 file
        return:
        -------
            refseq_genename_geneobject: structure:
                [gene_name] = list_of_geneobject_associated_with genename
            refseq_genename_geneid_transobject: structure:
                [gene_name][geneid][transid]= list_of_transobject_associated
                with the gene name, the geneid and the transid
            refseq_genename_geneid_transid_exonobject: structure:
                [gene_name][geneid][transid] = list_of_exons_objects
            refseq_genename_geneid_transid_cdsobject: structure:
                [genename][geneid][transid] = list of cds object
                associated with the transcript
            gene_name_to_synonyms_dict: dictionary that given a gene symbol,
                it returns the synonyms associated:
                gene_name_to_synonyms_dict[gene_name] = set of synonyms
    '''
    logging.info("Creating RefseqTranscript dataclasses")

    refseq_genename_geneobject = dict()
    refseq_genename_geneid_transobject = dict()
    refseq_genename_geneid_transid_exonobject = dict()
    refseq_genename_geneid_transid_cdsobject = dict()
    gene_name_to_synonym_dict = dict()
    current_gene_obj = None
    current_trans_obj = None
    current_exon_obj = None
    current_cds_obj = None
    gene_name = None

    genes_types = {
        'repeat_region',
        "enhancer",
        "silencer",
        "biological_region",
        'tandem_repeat',
        'origin_of_replication',
        'transcriptional_cis_regulatory_region',
        'sequence_comparison',
        'imprinting_control_region',
        'chromosome_breakpoint',
        'dnasei_hypersensitive_site',
        'insulator',
        'caat_signal',
        'epigenetically_modified_region',
        'd_loop',
        'protein_binding_site',
        'minisatellite',
        'meiotic_recombination_region',
        'nucleotide_cleavage_site',
        'replication_start_site',
        'regulatory_region',
        'nucleotide_motif',
        'tata_box',
        'conserved_region',
        'non_allelic_homologous_recombination_region',
        'microsatellite',
        'gc_rich_promoter_region',
        'gene',
        'response_element',
        'mobile_genetic_element',
        'mitotic_recombination_region',
        'sequence_secondary_structure',
        'matrix_attachment_site',
        'replication_regulatory_region',
        'sequence_alteration',
        'dispersed_repeat',
        'locus_control_region',
        'cage_cluster',
        'direct_repeat',
        'repeat_instability_region',
        'enhancer_blocking_element',
        'promoter',
        'pseudogene'
    }
    # Parsing transcripts into RefseqTranscripts dataclass
    transcripts_types = [
        "transcript",
        "rrna",
        "ncrna",
        "primary_transcript",
        "antisense_rna",
        "lnc_rna",
        "mirna",
        "snorna",
        "mrna",
        "snrna",
        "trna",
        "snoRNA"
    ]


    with gzip.open(refseq_gff, 'rt') as fin:
        for line in fin:
            gene_synonyms = list()
            feature_dict = dict()
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            feature = tmp[2].lower()
            refseq_chr = tmp[0]
            chr = convert_chromosomes(refseq_chr)
            pos = tmp[3]
            end = tmp[4]
            info = tmp[8].split(";")
            feature_dict["chromosome"] = chr
            feature_dict["start"] = pos
            feature_dict["end"] = end

            # parsing gene info
            if feature in genes_types:
                gene_id = None
                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    if field == "dbxref":
                        dbxref_results = result.split(",")
                        for dbxref_result in dbxref_results:
                            dbxref_result = dbxref_result.lower()
                            if "geneid" in dbxref_result:
                                gene_id = dbxref_result.split(":")[1]
                                feature_dict["gene_id"] = gene_id
                            elif "hgnc" in dbxref_result:
                                hgnc_id = dbxref_result.split(":")[-1]
                                hgnc_id = f"HGNC:{hgnc_id}"
                                feature_dict["hgnc_id"] = hgnc_id
                    elif field == "name":
                        gene_name = result
                    elif field == "gene_synonym":
                        if "," in result:
                            gene_synonyms = result.split(",")
                        else:
                            gene_synonyms.append(result)
                    feature_dict[field] = result

                if "name" not in feature_dict:
                    feature_dict["name"] = feature

                if "gene_biotype" not in feature_dict:
                    feature_dict["gene_biotype"] = feature

                # Adding gene synonyms to the dict
                if gene_synonyms:
                    for gene_synonym in gene_synonyms:
                        # creating a dict that maps gene_name : gene_synonyms
                        if gene_name in gene_name_to_synonym_dict:
                            gene_name_to_synonym_dict[gene_name].add(gene_synonym)
                        else:
                            gene_name_to_synonym_dict[gene_name] = {gene_synonym, gene_name}
                        # creating a dict that maps gene_synonym : gene_names
                        if gene_synonym in gene_name_to_synonym_dict:
                            gene_name_to_synonym_dict[gene_synonym].add(gene_name)
                            gene_name_to_synonym_dict[gene_synonym].update(gene_synonyms)
                        else:
                            gene_name_to_synonym_dict[gene_synonym] = {gene_name}
                            gene_name_to_synonym_dict[gene_synonym].update(gene_synonyms)

                # There is one gene that doesn't contain gene_id
                if "gene_id" not in feature_dict:
                    continue
                
                # just include the key:value pairs in feature dict that are dfined in the dataclass RefseqGene
                current_gene_obj = RefseqGene(**{k: v for k, v in feature_dict.items() if k in RefseqGene.__annotations__})

                if current_gene_obj is not None:
                    if gene_name not in refseq_genename_geneobject:
                        refseq_genename_geneobject[gene_name] = list()
                        refseq_genename_geneobject[gene_name].append(current_gene_obj)
                    else:
                        refseq_genename_geneobject[gene_name].append(current_gene_obj)
            
            elif feature in transcripts_types:
                current_trans_obj = None
                # set at the beggining because we have to increment the number of exons and cds 
                # for each transcript and the object 
                if current_gene_obj is not None:
                    current_gene_obj.increment_num_transcripts()
                
                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    if field == "parent":
                        if "-" in result:
                        # removing gene- from the gene name
                            gene_name = result.split("-", 1)[-1]
                    elif field == "gbkey":
                        feature_dict["transcript_type"] = result
                    elif field == "gene":
                        feature_dict["gene_name"] = result
                    elif field == "id":
                        # in V_gene_segment feature, it does not contain rna
                        if "rna-" not in result:
                            result.replace("id", "")
                        result = result.replace("rna-", "")
                        id = result
                        feature_dict["id"] = result
                    elif field == "transcript_id":
                        transcript_id = result
                        feature_dict["transcript_id"] = transcript_id
                        if transcript_id in mane_clin_refseq_id:
                            feature_dict["mane_clin"] = True
                        if transcript_id in mane_select_refseq_id:
                            feature_dict["mane_select"] = True
                        if transcript_id in lrg_refseq:
                            feature_dict["lrg_id"] = lrg_refseq[transcript_id]["lrg_id"]
                            feature_dict["lrg_transcript"] = lrg_refseq[transcript_id]["lrg_transcript"]
                            feature_dict["ccds"] = lrg_refseq[transcript_id]["ccds"]
                    elif field == "dbxref":
                        dbxref_results = result.split(",")
                        for dbxref_result in dbxref_results:
                            if "GeneID" in dbxref_result:
                                gene_id = dbxref_result.split(":")[1]
                                feature_dict["gene_id"] = gene_id
                            elif "hgnc" in dbxref_result:
                                hgnc_id = dbxref_result.split(":")[-1]
                                hgnc_id = f"HGNC:{hgnc_id}"
                                feature_dict["hgnc_id"] = hgnc_id
                if "name" not in feature_dict:
                    feature_dict["name"] = feature
                if "gene_biotype" not in feature_dict:
                    feature_dict["gene_biotype"] = feature
                
                # miRNA does not contain transcripts_ids so we will use it's id to identify them 
                if "transcript_id" not in feature_dict:
                    feature_dict["transcript_id"] = id
                
                # There are 13 transcript that don't contain transcript_id, to solve it:

                # if "transcript_id" not in feature_dict:
                #     continue
                current_trans_obj = RefseqTranscript(**{k: v for k, v in feature_dict.items() if k in RefseqTranscript.__annotations__})
                if current_trans_obj is not None:
                    if gene_name not in refseq_genename_geneid_transobject:
                        refseq_genename_geneid_transobject.setdefault(gene_name, dict())
                    if gene_id not in refseq_genename_geneid_transobject[gene_name]:
                        refseq_genename_geneid_transobject[gene_name].setdefault(gene_id, list())
                        refseq_genename_geneid_transobject[gene_name][gene_id].append(current_trans_obj)
                    else:
                        refseq_genename_geneid_transobject[gene_name][gene_id].append(current_trans_obj)

            elif feature == "exon":
                current_exon_obj = None
                if current_trans_obj is None:
                    continue
                if current_trans_obj is not None:
                    current_trans_obj.increment_numb_exons()
                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    if field == "id":
                        exon_id = result.replace("exon-", "")
                        feature_dict["exon_id"] = exon_id
                        if "-" in exon_id:
                            feature_dict["exon_number"] = exon_id.split("-")[1]
                        else:
                            feature_dict["exon_number"] = None

                    elif field == "gene":
                        gene_name = result
                    elif field == "parent":
                        transcript_id = result
                        # if there is no transcript_id field, we will use the parent as transcript_id
                        feature_dict["transcript_id"] = transcript_id
                    elif field == "dbxref":
                        dbxref_results = result.split(",")
                        for dbxref_result in dbxref_results:
                            if "GeneID" in dbxref_result:
                                gene_id = dbxref_result.split(":")[1]
                                feature_dict["id"] = gene_id
                    elif field == "gene_biotype":
                        field = "gene_type"
                    elif field == "gene":
                        id = result
                    feature_dict[field] = result
                # Adding only feature dict key:values which keys are defined
                # in RefseqExons dataclass
                current_exon_obj = RefseqExons(**{k: v for k, v in feature_dict.items() if k in RefseqExons.__annotations__})
                if current_exon_obj is not None:
                    if gene_name not in refseq_genename_geneid_transid_exonobject:
                        refseq_genename_geneid_transid_exonobject.setdefault(gene_name, dict())
                    if gene_id not in refseq_genename_geneid_transid_exonobject[gene_name]:
                        refseq_genename_geneid_transid_exonobject[gene_name].setdefault(gene_id, dict())
                    if current_trans_obj not in refseq_genename_geneid_transid_exonobject[gene_name][gene_id]:
                        refseq_genename_geneid_transid_exonobject[gene_name][gene_id].setdefault(current_trans_obj, list())
                        refseq_genename_geneid_transid_exonobject[gene_name][gene_id][current_trans_obj].append(current_exon_obj)
                    else:
                        refseq_genename_geneid_transid_exonobject[gene_name][gene_id][current_trans_obj].append(current_exon_obj)

            elif feature == "cds":
                current_cds_obj = None
                if current_trans_obj is None:
                    continue
                if current_trans_obj is not None:
                    current_trans_obj.increment_numb_cds()
                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    if field == "id":
                        feature_dict["id"] = result
                    elif field == "gene":
                        gene_name = result
                    elif field == "parent":
                        transcript_id = result.split("-")[1]
                        feature_dict["transcript_id"] = transcript_id
                    elif field == "dbxref":
                        dbxref_results = result.split(",")
                        for dbxref_result in dbxref_results:
                            if "GeneID" in dbxref_result:
                                gene_id = dbxref_result.split(":")[1]
                                feature_dict["id"] = gene_id
                feature_dict[field] = result
                
                current_cds_obj = RefseqCds(**{k: v for k, v in feature_dict.items() if k in RefseqCds.__annotations__})
                if current_cds_obj is not None:
                    if gene_name not in refseq_genename_geneid_transid_cdsobject:
                        refseq_genename_geneid_transid_cdsobject.setdefault(gene_name, dict())
                    if gene_id not in refseq_genename_geneid_transid_cdsobject[gene_name]:
                        refseq_genename_geneid_transid_cdsobject[gene_name].setdefault(gene_id, dict())
                    if transcript_id not in refseq_genename_geneid_transid_cdsobject[gene_name][gene_id]:
                        refseq_genename_geneid_transid_cdsobject[gene_name][gene_id].setdefault(current_trans_obj, list())
                        refseq_genename_geneid_transid_cdsobject[gene_name][gene_id][current_trans_obj].append(current_cds_obj)
                    else:
                        refseq_genename_geneid_transid_cdsobject[gene_name][gene_id][current_trans_obj].append(current_cds_obj)

            # if the element type is not gene_types, transcript_types, exon or cds
            # all current objects will be set to None. This avoids counting exons or cds
            # of a transcripts that does not correspond to
            else:
                current_cds_obj = None
                current_exon_obj = None
                current_gene_obj = None
                current_trans_obj = None

    return (
        refseq_genename_geneobject,
        refseq_genename_geneid_transobject,
        refseq_genename_geneid_transid_exonobject,
        refseq_genename_geneid_transid_cdsobject,
        gene_name_to_synonym_dict
    )

(
    refseq_genename_geneobject,
    refseq_genename_geneid_transobject,
    refseq_genename_geneid_transid_exonobject,
    refseq_genename_geneid_transid_cdsobject,
    gene_name_to_synonym_dict
) = get_refseq_transcripts(refseq_gff, mane_clin_refseq_id, mane_select_refseq_id)


