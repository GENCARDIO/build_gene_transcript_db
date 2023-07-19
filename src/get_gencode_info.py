from dataclasses import dataclass, field
import gzip
from collections import defaultdict
import logging

from src.get_mane import mane_clin_ens_id, mane_select_ens_id
from src.global_variables import logging
from src.get_lrg_info import lrg_ensembl

gencode_gff = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/gencode/gencode.v43.annotation.gff3.gz"

@dataclass
class EnsemblGenes:
    gene_id:                    str
    start:                      int
    end:                        int
    chromosome:                 str
    id:                         str
    gene_name:                  str
    gene_type:                  str
    level:                      int
    gene_synonyms:              str = field(default=None)
    tag:                        str = field(default=None)
    havana_gene:                str = field(default=None)
    hgnc_id:                    int = field(default=None)
    num_gencode_transcripts:    int = field(default=0)
    num_refseq_transcripts:     int = field(default=0)
    refseq_same_gene_coords:    str = field(default=None)

    def increment_num_gencode_transcripts(self):
        self.num_gencode_transcripts += 1
    
    def increment_num_refseq_transcripts(self):
        self.num_refseq_transcripts += 1
    
    def get_coords(self):
        return f"{self.chromosome}:{self.start}-{self.end}"


@dataclass
class EnsemblTranscript:
    transcript_id:              str
    start:                      int
    end:                        int
    chromosome:                 str
    transcript_type:            str
    transcript_name:            str
    id:                         str
    gene_id:                    str
    protein_id:                 str = field(default=None)
    lrg_transcript:             str = field(default=None)
    lrg_id:                     str = field(default=None)
    ccds:                       str = field(default=None)
    numb_exons:                 int = field(default=0)
    numb_cds:                   int = field(default=0)
    mane_select:                bool = field(default=False)
    mane_clin:                  bool = field(default=False)
    havana_transcript:          str = field(default=None)
    refseq_same_trans_coords:   str = field(default=None)
    refseq_same_exons_coords:   str = field(default=None)
    refseq_same_cds_coords:     str = field(default=None)

    def increment_numb_exons(self):
        self.numb_exons += 1

    def increment_numb_cds(self):
        self.numb_cds += 1

    def get_coords(self):
        return f"{self.chromosome}:{self.start}-{self.end}"

    def __hash__(self):
        return hash(self.transcript_id)

    def __eq__(self, other):
        if isinstance(other, EnsemblTranscript):
            return self.transcript_id == other.transcript_id
        return False


@dataclass
class EnsemblExons:
    exon_id:        str
    id:             str
    exon_number:    int
    start:          int
    end:            int
    chromosome:     str
    transcript_id:  str

    def get_coords(self):
        return f"{self.chromosome}:{self.start}-{self.end}"


@dataclass
class EnsemblCds:
    id:             str
    start:          int
    end:            int
    chromosome:     str
    exon_id:        str
    transcript_id:  str
    exon_number:    int

    def get_coords(self):
        return f"{self.chromosome}:{self.start}-{self.end}"


def parse_gene_gff3_line(feature_dict, gencode_genename_geneobject):

    for inf in feature_dict["info"]:
        field = inf.split("=")[0].lower()
        result = inf.split("=")[1]

        feature_dict[field] = result
        gene_name = feature_dict["gene_name"]

    # creting dataclass instance of EnsemblGenes with information of gene_dict that is
    # which its keys are defined as variables in the dataclass
    current_ensembl_gene = EnsemblGenes(**{k: v for k, v in gene_dict.items() if k in EnsemblGenes.__annotations__})
    if current_ensembl_gene is not None:
        if gene_name not in gencode_genename_geneobject:
            gencode_genename_geneobject[gene_name] = list()
            gencode_genename_geneobject[gene_name].append(current_ensembl_gene)
        else:
            gencode_genename_geneobject[gene_name].append(current_ensembl_gene)


def parse_gencode(
        gencode_gff: str,
        mane_clin_ens_id: set,
        mane_select_ens_id: set,
        lrg_ensembl: dict
):
    logging.info(f"Extracting information from {gencode_gff}")

    current_ensembl_gene = None
    current_ensembl_transcript = None
    current_ensembl_exon = None
    current_ensembl_cds = None

    gencode_genename_geneobject = dict()
    gencode_genename_geneid_transobject = dict()
    
    # Creates a dictionary with a default value of an empty dictionary
    gencode_genename_geneid_transid_exonobject = defaultdict(dict)
    gencode_genename_geneid_transid_cdsobject = defaultdict(dict)

    with gzip.open(gencode_gff, 'rt') as fin:
        for line in fin:
            feature_dict = dict()
            gene_dict = {}
            line = line.rstrip("\n")

            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            feature_dict["feature"] = tmp[2].lower()
            feature_dict["chromosome"] = tmp[0].replace("chr", "")
            feature_dict["start"] = tmp[3]
            feature_dict["end"] = tmp[4]
            feature_dict["info"] = tmp[8].split(";")

            # parsing genes and creating a EnsemblGene dataclass with its information
            if feature_dict["feature"] == "gene":
                for inf in feature_dict["info"]:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]

                    if field == "gene_name":
                        gene_name = result

                    elif field == "gene_id":
                        gene_id = result

                    feature_dict[field] = result

                # creting dataclass instance of EnsemblGenes with information of gene_dict that is
                # which its keys are defined as variables in the dataclass
                current_ensembl_gene = EnsemblGenes(**{k: v for k, v in gene_dict.items() if k in EnsemblGenes.__annotations__})
                if current_ensembl_gene is not None:
                    if gene_name not in gencode_genename_geneobject:
                        gencode_genename_geneobject[gene_name] = list()
                        gencode_genename_geneobject[gene_name].append(current_ensembl_gene)
                    else:
                        gencode_genename_geneobject[gene_name].append(current_ensembl_gene)

            # parsing transcripts and creating a EnsemblTranscript dataclass with its information
            elif feature == "transcript":

                # Incrementing the attribute num_gencode_transcripts by 1 (numb of transcripts that the current gene have + 1)
                current_ensembl_gene.increment_num_gencode_transcripts()
                transcript_dict = dict()
                transcript_dict["chromosome"] = chr
                transcript_dict["start"] = pos
                transcript_dict["end"] = end

                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    feature_fields.add(field)

                    if field == "transcript_id":
                        transcript_id = result
                        if transcript_id in lrg_ensembl:
                            transcript_dict["lrg_id"] = lrg_ensembl[transcript_id]["lrg_id"]
                            transcript_dict["lrg_transcript"] = lrg_ensembl[transcript_id]["lrg_transcript"]
                            transcript_dict["ccds"] = lrg_ensembl[transcript_id]["ccds"]

                    elif field == "gene_name":
                        gene_name = result
                    
                    elif field == "gene_id":
                        gene_id = result
                    
                    # adding all transcript information in transcript_dict
                    transcript_dict[field] = result
                    # adding mane_select and clin
                if transcript_id in mane_select_ens_id:
                    transcript_dict["mane_select"] = True
                if transcript_id in mane_clin_ens_id:
                    transcript_dict["mane_clin"] = True

                current_ensembl_transcript = EnsemblTranscript(**{k: v for k, v in transcript_dict.items() if k in EnsemblTranscript.__annotations__})
                if current_ensembl_transcript is not None:
                    if gene_name not in gencode_genename_geneid_transobject:
                        gencode_genename_geneid_transobject.setdefault(gene_name, dict())
                    if gene_id not in  gencode_genename_geneid_transobject[gene_name]:
                        gencode_genename_geneid_transobject[gene_name].setdefault(gene_id, list())
                        gencode_genename_geneid_transobject[gene_name][gene_id].append(current_ensembl_transcript)
                    else:
                        gencode_genename_geneid_transobject[gene_name][gene_id].append(current_ensembl_transcript)

            # parsing exons and creating a EnsemblExon dataclass with its information
            elif feature == "exon":
                current_ensembl_transcript.increment_numb_exons()

                exon_dict = dict()
                exon_dict["chromosome"] = chr
                exon_dict["start"] = pos
                exon_dict["end"] = end

                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    feature_fields.add(field)

                    if field == "gene_name":
                        gene_name = result
                    
                    elif field == "gene_id":
                        gene_id = result

                    exon_dict[field] = result

                current_ensembl_exon = EnsemblExons(**{k: v for k, v in exon_dict.items() if k in EnsemblExons.__annotations__})
                if current_ensembl_exon is not None:
                    if gene_name not in gencode_genename_geneid_transid_exonobject:
                        gencode_genename_geneid_transid_exonobject.setdefault(gene_name, dict())
                    if gene_id not in  gencode_genename_geneid_transid_exonobject[gene_name]:
                        gencode_genename_geneid_transid_exonobject[gene_name].setdefault(gene_id, dict())
                    if current_ensembl_transcript not in gencode_genename_geneid_transid_exonobject[gene_name][gene_id]:
                        gencode_genename_geneid_transid_exonobject[gene_name][gene_id].setdefault(current_ensembl_transcript, list())
                        gencode_genename_geneid_transid_exonobject[gene_name][gene_id][current_ensembl_transcript].append(current_ensembl_exon)
                    else:
                        gencode_genename_geneid_transid_exonobject[gene_name][gene_id][current_ensembl_transcript].append(current_ensembl_exon)
            
            # parsing cds and creating a EnsemblCds dataclass with its information
            elif feature == "cds":
                current_ensembl_transcript.increment_numb_cds()

                cds_dict = dict()
                cds_dict["chromosome"] = chr
                cds_dict["start"] = pos
                cds_dict["end"] = end

                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    feature_fields.add(field)

                    if field == "id":
                        cds_id = result

                    elif field == "gene_name":
                        gene_name = result
                    
                    elif field == "parent":
                        transcript_id = result
                    
                    elif field == "gene_id":
                        gene_id = result

                    cds_dict[field] = result

                current_ensembl_cds = EnsemblCds(**{k: v for k, v in cds_dict.items() if k in EnsemblCds.__annotations__})

                if current_ensembl_cds is not None:
                    if gene_name not in gencode_genename_geneid_transid_cdsobject:
                        gencode_genename_geneid_transid_cdsobject.setdefault(gene_name, dict())
                    if gene_id not in gencode_genename_geneid_transid_cdsobject[gene_name]:
                        gencode_genename_geneid_transid_cdsobject[gene_name].setdefault(gene_id, dict())
                    if transcript_id not in gencode_genename_geneid_transid_cdsobject[gene_name][gene_id]:
                        gencode_genename_geneid_transid_cdsobject[gene_name][gene_id].setdefault(current_ensembl_transcript, list())
                        gencode_genename_geneid_transid_cdsobject[gene_name][gene_id][current_ensembl_transcript].append(current_ensembl_cds)
                    else:
                        gencode_genename_geneid_transid_cdsobject[gene_name][gene_id][current_ensembl_transcript].append(current_ensembl_cds)

    logging.info(f"{gencode_gff} has been parsed successfully")
    
    # Adding the last gene, transcript, exon and cds found in the gff file
    # as it is only added when the following line with the same element is found
    # as consequence the lasts components are not included in the dictionary until now

    return (gencode_genename_geneobject, gencode_genename_geneid_transobject, gencode_genename_geneid_transid_exonobject, gencode_genename_geneid_transid_cdsobject)


(
    gencode_genename_geneobject,
    gencode_genename_geneid_transobject,
    gencode_genename_geneid_transid_exonobject,
    gencode_genename_geneid_transid_cdsobject
) = parse_gencode(
    gencode_gff,
    mane_clin_ens_id,
    mane_select_ens_id,
    lrg_ensembl
)

# gencode_gene_file = "gencode_genename_geneobject.txt"
# gencode_transcript_file = "gencode_genename_geneid_transobject.txt"
# gencode_exon_file = "gencode_genename_geneid_transid_exonobject.txt"
# gencode_cds_file = "gencode_genename_geneid_transid_cdsobject.txt"

# with open(gencode_gene_file, "w") as f:
#     json.dump(gencode_genename_geneobject, f)
# logging.info("gencode_genename_geneobject.txt created successfully")

# with open(gencode_transcript_file, "w") as f2:
#     json.dump(gencode_genename_geneid_transobject, f2)
# logging.info("gencode_genename_geneid_transobject.txt created successfully")

# with open(gencode_exon_file, "w") as f3:
#     json.dump(gencode_genename_geneid_transid_exonobject, f3)
# logging.info("gencode_genename_geneid_transid_exonobject.txt created successfully")

# with open(gencode_cds_file, "w") as f4:
#     json.dump(gencode_genename_geneid_transid_cdsobject, f4)
# logging.info("gencode_genename_geneid_transid_cdsobject.txt created successfully")
