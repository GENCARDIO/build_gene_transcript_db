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


def parse_gff3_line_general_info(line):
    feature_dict = dict()
    line = line.rstrip("\n")

    if line.startswith("#"):
        return None
    tmp = line.split("\t")
    feature_dict["feature"] = tmp[2].lower()
    feature_dict["chromosome"] = tmp[0].replace("chr", "")
    feature_dict["start"] = tmp[3]
    feature_dict["end"] = tmp[4]
    feature_dict["info"] = tmp[8].split(";")

    return (feature_dict)


def parse_gene_gff3_line(
    feature_dict: dict,
    gencode_genename_geneobject: dict
):

    for inf in feature_dict["info"]:
        field = inf.split("=")[0].lower()
        result = inf.split("=")[1]
        feature_dict[field] = result
    gene_name = feature_dict["gene_name"]

    # creting dataclass instance of EnsemblGenes with information of gene_dict that is
    # which its keys are defined as variables in the dataclass
    current_ensembl_gene = EnsemblGenes(**{k: v for k, v in feature_dict.items() if k in EnsemblGenes.__annotations__})
    if current_ensembl_gene is not None:
        if gene_name not in gencode_genename_geneobject:
            gencode_genename_geneobject[gene_name] = list()
            gencode_genename_geneobject[gene_name].append(current_ensembl_gene)
        else:
            gencode_genename_geneobject[gene_name].append(current_ensembl_gene)

    return (current_ensembl_gene, gencode_genename_geneobject)


def parse_transcript_gff3_line(
    feature_dict: dict,
    gencode_genename_geneid_transobject,
    mane_select_ens_id,
    mane_clin_ens_id,
    lrg_ensembl
):
    for inf in feature_dict["info"]:
        field = inf.split("=")[0].lower()
        result = inf.split("=")[1]

        if field == "transcript_id":
            transcript_id = result
            if transcript_id in lrg_ensembl:
                feature_dict["lrg_id"] = lrg_ensembl[transcript_id]["lrg_id"]
                feature_dict["lrg_transcript"] = lrg_ensembl[transcript_id]["lrg_transcript"]
                feature_dict["ccds"] = lrg_ensembl[transcript_id]["ccds"]

        elif field == "gene_name":
            gene_name = result

        elif field == "gene_id":
            gene_id = result

        # adding all transcript information in feature_dict
        feature_dict[field] = result
        # adding mane_select and clin
    if transcript_id in mane_select_ens_id:
        feature_dict["mane_select"] = True
    if transcript_id in mane_clin_ens_id:
        feature_dict["mane_clin"] = True

    current_ensembl_transcript = EnsemblTranscript(**{k: v for k, v in feature_dict.items() if k in EnsemblTranscript.__annotations__})
    if current_ensembl_transcript is not None:
        if gene_name not in gencode_genename_geneid_transobject:
            gencode_genename_geneid_transobject.setdefault(gene_name, dict())
        if gene_id not in gencode_genename_geneid_transobject[gene_name]:
            gencode_genename_geneid_transobject[gene_name].setdefault(gene_id, list())
            gencode_genename_geneid_transobject[gene_name][gene_id].append(current_ensembl_transcript)
        else:
            gencode_genename_geneid_transobject[gene_name][gene_id].append(current_ensembl_transcript)

    return (current_ensembl_transcript, gencode_genename_geneid_transobject)


def parse_exon_gff3_line(
    feature_dict,
    current_ensembl_transcript,
    gencode_genename_geneid_transid_exonobject
):

    for inf in feature_dict["info"]:
        field = inf.split("=")[0].lower()
        result = inf.split("=")[1]

        if field == "gene_name":
            gene_name = result
        
        elif field == "gene_id":
            gene_id = result

        feature_dict[field] = result

    current_ensembl_exon = EnsemblExons(**{k: v for k, v in feature_dict.items() if k in EnsemblExons.__annotations__})
    if current_ensembl_exon is not None:
        if gene_name not in gencode_genename_geneid_transid_exonobject:
            gencode_genename_geneid_transid_exonobject.setdefault(gene_name, dict())
        if gene_id not in gencode_genename_geneid_transid_exonobject[gene_name]:
            gencode_genename_geneid_transid_exonobject[gene_name].setdefault(gene_id, dict())
        if current_ensembl_transcript not in gencode_genename_geneid_transid_exonobject[gene_name][gene_id]:
            gencode_genename_geneid_transid_exonobject[gene_name][gene_id].setdefault(current_ensembl_transcript, list())
            gencode_genename_geneid_transid_exonobject[gene_name][gene_id][current_ensembl_transcript].append(current_ensembl_exon)
        else:
            gencode_genename_geneid_transid_exonobject[gene_name][gene_id][current_ensembl_transcript].append(current_ensembl_exon)

    return (gencode_genename_geneid_transid_exonobject)


def parse_cds_gff3_line(
    feature_dict,
    current_ensembl_transcript,
    gencode_genename_geneid_transid_cdsobject
):
    for inf in feature_dict["info"]:
        field = inf.split("=")[0].lower()
        result = inf.split("=")[1]

        if field == "gene_name":
            gene_name = result

        elif field == "parent":
            transcript_id = result

        elif field == "gene_id":
            gene_id = result

        feature_dict[field] = result

    current_ensembl_cds = EnsemblCds(**{k: v for k, v in feature_dict.items() if k in EnsemblCds.__annotations__})

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

    return (gencode_genename_geneid_transid_cdsobject)


def parse_gencode(
        gencode_gff: str,
        mane_clin_ens_id: set,
        mane_select_ens_id: set,
        lrg_ensembl: dict
):
    logging.info(f"Extracting information from {gencode_gff}")

    current_ensembl_gene = None
    current_ensembl_transcript = None

    gencode_genename_geneobject = dict()
    gencode_genename_geneid_transobject = dict()

    # Creates a dictionary with a default value of an empty dictionary
    gencode_genename_geneid_transid_exonobject = defaultdict(dict)
    gencode_genename_geneid_transid_cdsobject = defaultdict(dict)

    with gzip.open(gencode_gff, 'rt') as fin:
        for line in fin:
            feature_dict = parse_gff3_line_general_info(line)
            # when line starts with # it is returned None and we
            # won't analyse the line
            if feature_dict is None:
                continue

            # parse genes and create EnsemblGene dataclass with its information
            if feature_dict["feature"] == "gene":
                (
                    current_ensembl_gene,
                    gencode_genename_geneobject
                ) = parse_gene_gff3_line(
                    feature_dict,
                    gencode_genename_geneobject
                )
            # parse transcript and create EnsemblTranscript instance
            # with its information
            elif feature_dict["feature"] == "transcript":

                # Incrementing the attribute num_gencode_transcripts by 1
                # (numb of transcripts that the current gene have + 1)
                current_ensembl_gene.increment_num_gencode_transcripts()

                (
                    current_ensembl_transcript,
                    gencode_genename_geneid_transobject
                ) = parse_transcript_gff3_line(
                    feature_dict,
                    gencode_genename_geneid_transobject,
                    mane_select_ens_id,
                    mane_clin_ens_id,
                    lrg_ensembl
                )

            # parse exon and create EnsemblExon dataclass with its information
            elif feature_dict["feature"] == "exon":
                current_ensembl_transcript.increment_numb_exons()
                gencode_genename_geneid_transid_exonobject = (
                    parse_exon_gff3_line(
                        feature_dict,
                        current_ensembl_transcript,
                        gencode_genename_geneid_transid_exonobject
                    )
                )

            # parse cds and create EnsemblCds dataclass with its information
            elif feature_dict["feature"] == "cds":
                current_ensembl_transcript.increment_numb_cds()
                gencode_genename_geneid_transid_cdsobject = (
                    parse_cds_gff3_line(
                        feature_dict,
                        current_ensembl_transcript,
                        gencode_genename_geneid_transid_cdsobject
                    )
                )

    logging.info(f"{gencode_gff} has been parsed successfully")

    return (
        gencode_genename_geneobject,
        gencode_genename_geneid_transobject,
        gencode_genename_geneid_transid_exonobject,
        gencode_genename_geneid_transid_cdsobject
    )


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
