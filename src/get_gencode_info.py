from dataclasses import dataclass, field
import gzip
from collections import defaultdict
import logging
from typing import ClassVar

from src.get_mane import get_mane, mane_gff
from src.global_variables import (
    logging,
    lo_38_to_37,  # To perform liftover from grch38 to 37
    lo_37_to_38  # To perform liftover from grch37 to 38
) 
from src.get_lrg_info import get_lrg_trancripts, lrg_gff

gencode_gff_grch38 = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/gencode/gencode.v43.annotation.gff3.gz"
gencode_gff_grch37 = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/gencode/gencode.v44lift37.annotation.gff3.gz"
mock_gencode_gff_grch38 = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/gencode/mock.gencode.v43.annotation.gff3.gz"
mock_gencode_gff_grch37 = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/gencode/mock.gencode.v44lift37.annotation.gff3.gz"

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
    genome_version:             int
    liftover_chromosome:        str = field(default=None)
    liftover_start:             int = field(default=None)
    liftover_end:               int = field(default=None)
    gene_synonyms:              str = field(default=None)
    tag:                        str = field(default=None)
    havana_gene:                str = field(default=None)
    hgnc_id:                    int = field(default=None)
    num_gencode_transcripts:    int = field(default=0)
    refseq_same_gene_coords:    str = field(default=None)
    transcripts:                str = field(default=None)

    def get_full_id(self):
        return (f"{self.id}_{self.genome_version}")
    
    def __hash__(self):
        return hash(self.gene_id)
    
    def __eq__(self, other):
        if isinstance(other, EnsemblGenes):
            return self.gene_id == other.gene_id
        return False

    def increment_num_gencode_transcripts(self):
        self.num_gencode_transcripts += 1
    
    def add_transcript(self, transcript_obj):
        # check if exon transcript_id and transcript object id is the same:
        if transcript_obj.gene_id == self.gene_id:
            if self.transcripts is None:
                self.transcript = [transcript_obj]
            else:
                self.transcript.append[transcript_obj]
        else:
            logging.critic(f"transcript: {transcript_obj} is not associated \
                with gene: {self}")

    def get_coords(self):
        return f"{self.chromosome}:{int(self.start)}-{int(self.end)}"

    def get_coordinate_liftover(self):
        chromosome_name = f"chr{self.chromosome}"

        if self.genome_version == 38:    
            start = lo_38_to_37.convert_coordinate(
                chromosome_name,
                int(self.start)
            )
            end = lo_38_to_37.convert_coordinate(
                chromosome_name,
                int(self.end)
            )
            if start:
                # liftover performed successfully
                (
                    self.liftover_chromosome,
                    self.liftover_start,
                    strand,
                    _
                ) = start[0]

            if end:
                (
                    _,
                    self.liftover_end,
                    strand,
                    _
                ) = end[0]
                logging.info(
                    f"liftover performed successfully for {self.gene_id}: \
                    {self.liftover_chromosome}:{self.liftover_start}-{self.liftover_end}"
                )
        elif self.genome_version == 37:

            start = lo_37_to_38.convert_coordinate(
                chromosome_name,
                int(self.start),
            )
            end = lo_37_to_38.convert_coordinate(
                chromosome_name,
                int(self.end)
            )
            if start:
                # liftover performed successfully
                (
                    self.liftover_chromosome,
                    self.liftover_start,
                    _,
                    _
                ) = start[0]
            if end:
                (
                    self.liftover_chromosome,
                    self.liftover_end,
                    _,
                    _
                ) = end[0]
                logging.info(
                    f"liftover performed successfully for {self.gene_id}\
                    {self.liftover_chromosome}:{self.liftover_start}-{self.liftover_end}"
                )


# @dataclass
# class EnsemblGenes_Grch37:
#     gene_id:                    str
#     start:                      int
#     end:                        int
#     chromosome:                 str
#     id:                         str
#     gene_name:                  str
#     gene_type:                  str
#     level:                      int
#     hg19_start:                 int = field(default=None)
#     hg19_end:                   int = field(default=None)
#     gene_synonyms:              str = field(default=None)
#     tag:                        str = field(default=None)
#     havana_gene:                str = field(default=None)
#     hgnc_id:                    int = field(default=None)
#     num_gencode_transcripts:    int = field(default=0)
#     refseq_same_gene_coords:    str = field(default=None)
#     transcripts:                str = field(default=None)
#     grch38_lo_chromosome:       int = field(default=None)
#     grch38_lo_start:            int = field(default=None)
#     grch38_lo_end:              int = field(default=None)

#     def __hash__(self):
#         return hash(self.gene_id)
    
#     def __eq__(self, other):
#         if isinstance(other, EnsemblGenes):
#             return self.gene_id == other.gene_id
#         return False

#     def increment_num_gencode_transcripts(self):
#         self.num_gencode_transcripts += 1
    
#     def add_transcript(self, transcript_obj):
#         # check if exon transcript_id and transcript object id is the same:
#         if transcript_obj.gene_id == self.gene_id:
#             if self.transcripts is None:
#                 self.transcript = [transcript_obj]
#             else:
#                 self.transcript.append[transcript_obj]
#         else:
#             logging.critic(f"transcript: {transcript_obj} is not associated \
#                 with gene: {self}")

#     def get_coords(self):
#         return f"{self.chromosome}:{int(self.start)}-{int(self.end)}"

#     def get_coordinate_liftover(self):
#         grch38 = lo_37_to_38.convert_coordinate(
#             self.chromosome,
#             int(self.start),
#             int(self.end)
#         )
#         if grch38:
#             # liftover performed successfully
#             (
#                 self.grch38_lo_chromosome,
#                 self.grch38_lo_start,
#                 self.grch38_lo_end,
#                 _
#             ) = grch38[0]
#             logging.info(f"liftover performed successfully for {self.gene_id}")


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
    genome_version:             int
    parent:                     str
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
    transcript_exons:           str = field(default=None)
    transcript_cds:             str = field(default=None)
    liftover_chromosome:       int = field(default=None)
    liftover_start:            int = field(default=None)
    liftover_end:              int = field(default=None)
    
    def get_full_id(self):
        return (f"{self.id}_{self.genome_version}")
    
    def increment_numb_exons(self):
        self.numb_exons += 1

    def increment_numb_cds(self):
        self.numb_cds += 1

    def add_exon(self, exon_obj):
        if exon_obj.transcript_id == self.transcript_id:
            if self.transcript_exons is None:
                self.transcript_exons = [exon_obj]
            else:
                self.transcript_exons.append(exon_obj)
        else:
            logging.critical(f"exon {exon_obj} not match transid with {self}")

    def add_cds(self, cds_obj):
        if cds_obj.transcript_id == self.transcript_id:
            if self.transcript_cds is None:
                self.transcript_cds = [cds_obj]
            else:
                self.transcript_cds.append(cds_obj)
        else:
            logging.critical(f"cds {cds_obj} not match transid with {self}")
            raise ValueError(f"cds {cds_obj.transcript_id} not match transid with {self.transcript_id}")

    def get_coords(self):
        return f"{self.chromosome}:{int(self.start)}-{int(self.end)}"

    def __hash__(self):
        return hash(self.transcript_id)

    def __eq__(self, other):
        if isinstance(other, EnsemblTranscript):
            return self.transcript_id == other.transcript_id
        return False

    def get_coordinate_liftover(self):
        chromosome_name = f"chr{self.chromosome}"

        if self.genome_version == 38:
            start = lo_38_to_37.convert_coordinate(
                chromosome_name,
                int(self.start)
            )
            end = lo_38_to_37.convert_coordinate(
                chromosome_name,
                int(self.end)
            )
            if start:
                # liftover performed successfully
                (
                    self.liftover_chromosome,
                    self.liftover_start,
                    strand,
                    _
                ) = start[0]

            if end:
                (
                    _,
                    self.liftover_end,
                    strand,
                    _
                ) = end[0]
                logging.info(
                    f"liftover performed successfully for {self.gene_id}\
                    {self.liftover_chromosome}:{self.liftover_start}-{self.liftover_end}"
                )
        elif self.genome_version == 37:
            start = lo_37_to_38.convert_coordinate(
                chromosome_name,
                int(self.start),
            )
            end = lo_37_to_38.convert_coordinate(
                chromosome_name,
                int(self.end)
            )
            if start:
                # liftover performed successfully
                (
                    self.liftover_chromosome,
                    self.liftover_start,
                    _,
                    _
                ) = start[0]
            if end:
                (
                    self.liftover_chromosome,
                    self.liftover_end,
                    _,
                    _
                ) = end[0]

                logging.info(
                    f"liftover performed successfully for {self.gene_id}\
                    {self.liftover_chromosome}:{self.liftover_start}-{self.liftover_end}"
                )


# @dataclass
# class EnsemblTranscript_Grch37:
#     transcript_id:              str
#     start:                      int
#     end:                        int
#     chromosome:                 str
#     transcript_type:            str
#     transcript_name:            str
#     id:                         str
#     gene_id:                    str
#     protein_id:                 str = field(default=None)
#     lrg_transcript:             str = field(default=None)
#     lrg_id:                     str = field(default=None)
#     ccds:                       str = field(default=None)
#     numb_exons:                 int = field(default=0)
#     numb_cds:                   int = field(default=0)
#     mane_select:                bool = field(default=False)
#     mane_clin:                  bool = field(default=False)
#     havana_transcript:          str = field(default=None)
#     refseq_same_trans_coords:   str = field(default=None)
#     refseq_same_exons_coords:   str = field(default=None)
#     refseq_same_cds_coords:     str = field(default=None)
#     transcript_exons:           str = field(default=None)
#     transcript_cds:             str = field(default=None)
#     grch38_lo_chromosome:       int = field(default=None)
#     grch38_lo_start:            int = field(default=None)
#     grch38_lo_end:              int = field(default=None)

#     def increment_numb_exons(self):
#         self.numb_exons += 1

#     def increment_numb_cds(self):
#         self.numb_cds += 1

#     def add_exon(self, exon_obj):
#         if exon_obj.transcript_id == self.transcript_id:
#             if self.transcript_exons is None:
#                 self.transcript_exons = [exon_obj]
#             else:
#                 self.transcript_exons.append(exon_obj)
#         else:
#             logging.critical(f"exon {exon_obj} not match transid with {self}")

#     def add_cds(self, cds_obj):
#         if cds_obj.transcript_id == self.transcript_id:
#             if self.transcript_cds is None:
#                 self.transcript_cds = [cds_obj]
#             else:
#                 self.transcript_cds.append(cds_obj)
#         else:
#             logging.critical(f"cds {cds_obj} not match transid with {self}")

#     def get_coords(self):
#         return f"{self.chromosome}:{int(self.start)}-{int(self.end)}"

#     def __hash__(self):
#         return hash(self.transcript_id)

#     def __eq__(self, other):
#         if isinstance(other, EnsemblTranscript):
#             return self.transcript_id == other.transcript_id
#         return False

#     def get_coordinate_liftover(self):
#         grch38 = lo_37_to_38.convert_coordinate(
#             self.chromosome,
#             int(self.start),
#             int(self.end)
#         )
#         if grch38:
#             # liftover performed successfully
#             (
#                 self.grch38_lo_chromosome,
#                 self.grch38_lo_start,
#                 self.grch38_lo_end,
#                 _
#             ) = grch38[0]
#             logging.info(f"liftover performed successfully for {self.gene_id}")

@dataclass
class EnsemblExons:
    exon_id:                str
    id:                     str
    exon_number:            int
    start:                  int
    end:                    int
    chromosome:             str
    transcript_id:          str
    genome_version:         int
    liftover_chromosome:    int = field(default=None)
    liftover_start:         int = field(default=None)
    liftover_end:           int = field(default=None)

    def get_full_id(self):
        return (f"{self.id}_{self.genome_version}")

    def get_coords(self):
        return f"{self.chromosome}:{int(self.start)}-{int(self.end)}"
    
    def get_coordinate_liftover(self):
        chromosome_name = f"chr{self.chromosome}"

        if self.genome_version == 38:
            start = lo_38_to_37.convert_coordinate(
                chromosome_name,
                int(self.start)
            )
            end = lo_38_to_37.convert_coordinate(
                chromosome_name,
                int(self.end)
            )
            if start:
                # liftover performed successfully
                (
                    self.liftover_chromosome,
                    self.liftover_start,
                    strand,
                    _
                ) = start[0]

            if end:
                (
                    _,
                    self.liftover_end,
                    strand,
                    _
                ) = end[0]
                logging.info(
                    f"liftover performed successfully for {self.exon_id}\
                    {self.liftover_chromosome}:{self.liftover_start}-{self.liftover_end}"
                )
        elif self.genome_version == 37:
            start = lo_37_to_38.convert_coordinate(
                chromosome_name,
                int(self.start),
            )
            end = lo_37_to_38.convert_coordinate(
                chromosome_name,
                int(self.end)
            )
            if start:
                # liftover performed successfully
                (
                    self.liftover_chromosome,
                    self.liftover_start,
                    _,
                    _
                ) = start[0]
            if end:
                (
                    self.liftover_chromosome,
                    self.liftover_end,
                    _,
                    _
                ) = end[0]

                logging.info(
                    f"liftover performed successfully for {self.exon_id}\
                    {self.liftover_chromosome}:{self.liftover_start}-{self.liftover_end}"
                    )


# @dataclass
# class EnsemblExons_Grch37:
#     exon_id:                str
#     id:                     str
#     exon_number:            int
#     start:                  int
#     end:                    int
#     chromosome:             str
#     transcript_id:          str
#     grch38_lo_chromosome:   int = field(default=None)
#     grch38_lo_start:        int = field(default=None)
#     grch38_lo_end:          int = field(default=None)
    
#     def get_coords(self):
#         return f"{self.chromosome}:{int(self.start)}-{int(self.end)}"

#     def get_coordinate_liftover(self):
#         grch38 = lo_37_to_38.convert_coordinate(
#             self.chromosome,
#             int(self.start),
#             int(self.end)
#         )
#         if grch38:
#             # liftover performed successfully
#             (
#                 self.grch38_lo_chromosome,
#                 self.grch38_lo_start,
#                 self.grch38_lo_end,
#                 _
#             ) = grch38[0]
#             logging.info(f"liftover performed successfully for {self.gene_id}")


@dataclass
class EnsemblCds:

    _last_primary_key: ClassVar[int] = 0
    id:                         str
    start:                      int
    end:                        int
    chromosome:                 str
    exon_id:                    str
    transcript_id:              str
    exon_number:                int
    genome_version:             int
    parent:                     str
    liftover_chromosome:        int = field(default=None)
    liftover_start:             int = field(default=None)
    liftover_end:               int = field(default=None)

    @classmethod
    def _get_next_primary_key(cls) -> int:
        cls._last_primary_key += 1
        return (cls._last_primary_key)
    
    def __post_init__(self):
        self.primary_key = self._get_next_primary_key()

    def get_full_id(self):
        return (f"{self.id}_{self.genome_version}")

    def get_coords(self):
        return f"{self.chromosome}:{int(self.start)}-{int(self.end)}"

    def get_coordinate_liftover(self):
        chromosome_name = f"chr{self.chromosome}"
        if self.genome_version == 38:
            start = lo_38_to_37.convert_coordinate(
                chromosome_name,
                int(self.start)
            )
            end = lo_38_to_37.convert_coordinate(
                chromosome_name,
                int(self.end)
            )
            if start:
                # liftover performed successfully
                (
                    self.liftover_chromosome,
                    self.liftover_start,
                    strand,
                    _
                ) = start[0]

            if end:
                (
                    _,
                    self.liftover_end,
                    strand,
                    _
                ) = end[0]
                logging.info(
                    f"liftover performed successfully for {self.id}\
                    {self.liftover_chromosome}:{self.liftover_start}-{self.liftover_end}"
                )
        elif self.genome_version == 37:
            start = lo_37_to_38.convert_coordinate(
                chromosome_name,
                int(self.start),
            )
            end = lo_37_to_38.convert_coordinate(
                chromosome_name,
                int(self.end)
            )
            if start:
                # liftover performed successfully
                (
                    self.liftover_chromosome,
                    self.liftover_start,
                    _,
                    _
                ) = start[0]
            if end:
                (
                    self.liftover_chromosome,
                    self.liftover_end,
                    _,
                    _
                ) = end[0]

                logging.info(
                    f"liftover performed successfully for {self.id}\
                    {self.liftover_chromosome}:{self.liftover_start}-{self.liftover_end}"
                )

# @dataclass
# class EnsemblCds_Grch37:
#     id:                     str
#     start:                  int
#     end:                    int
#     chromosome:             str
#     exon_id:                str
#     transcript_id:          str
#     exon_number:            int
#     grch38_lo_chromosome:   int = field(default=None)
#     grch38_lo_start:        int = field(default=None)
#     grch38_lo_end:          int = field(default=None)

#     def get_coords(self):
#         return f"{self.chromosome}:{int(self.start)}-{int(self.end)}"

#     def get_coordinate_liftover(self):
#         grch38 = lo_37_to_38.convert_coordinate(
#             self.chromosome,
#             int(self.start),
#             int(self.end)
#         )
#         if grch38:
#             # liftover performed successfully
#             (
#                 self.grch38_lo_chromosome,
#                 self.grch38_lo_start,
#                 self.grch38_lo_end,
#                 _
#             ) = grch38[0]
#             logging.info(f"liftover performed successfully for {self.gene_id}")



def parse_gff3_line_general_info(line):
    feature_dict = dict()
    line = line.rstrip("\n")

    if line.startswith("#"):
        return None
    tmp = line.split("\t")
    feature_dict["feature"] = tmp[2].lower()
    feature_dict["chromosome"] = str(tmp[0].replace("chr", ""))
    feature_dict["start"] = tmp[3]
    feature_dict["end"] = tmp[4]
    feature_dict["info"] = tmp[8].split(";")

    return (feature_dict)


def parse_gene_gff3_line(
    feature_dict: dict,
    gencode_genename_geneobject: dict,
    genome_v: int
):

    for inf in feature_dict["info"]:
        field = inf.split("=")[0].lower()
        result = inf.split("=")[1].upper()
        feature_dict[field] = result
    gene_name = feature_dict["gene_name"]
    
    if genome_v == 38:
        feature_dict["genome_version"] = 38
    elif genome_v == 37:
        feature_dict["genome_version"] = 37
    else:
        raise ValueError(f"genome version: {genome_v} is not 37 or 38")
    
    # creting dataclass instance of EnsemblGenes with information of gene_dict that is
    # which its keys are defined as variables in the dataclass
    current_ensembl_gene = EnsemblGenes(**{k: v for k, v in feature_dict.items() if k in EnsemblGenes.__annotations__})
    if current_ensembl_gene is not None:
        current_ensembl_gene.get_coordinate_liftover()
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
    lrg_ensembl,
    genome_v: int
):
    for inf in feature_dict["info"]:
        field = inf.split("=")[0].lower()
        result = inf.split("=")[1].upper()

        if field == "id":
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
    
    feature_dict["genome_version"] = genome_v

    current_ensembl_transcript = EnsemblTranscript(**{k: v for k, v in feature_dict.items() if k in EnsemblTranscript.__annotations__})

    if current_ensembl_transcript is not None:
        current_ensembl_transcript.get_coordinate_liftover()
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
    gencode_genename_geneid_transid_exonobject,
    genome_v
):

    for inf in feature_dict["info"]:
        field = inf.split("=")[0].lower()
        result = inf.split("=")[1].upper()

        if field == "gene_name":
            gene_name = result
        elif field == "id":
            result = result.replace("EXON:", "")
        elif field == "gene_id":
            gene_id = result

        feature_dict[field] = result
    feature_dict["genome_version"] = genome_v 
    current_ensembl_exon = EnsemblExons(**{k: v for k, v in feature_dict.items() if k in EnsemblExons.__annotations__})

    if current_ensembl_exon is not None:
        current_ensembl_exon.get_coordinate_liftover()
        if gene_name not in gencode_genename_geneid_transid_exonobject:
            gencode_genename_geneid_transid_exonobject.setdefault(gene_name, dict())
        if gene_id not in gencode_genename_geneid_transid_exonobject[gene_name]:
            gencode_genename_geneid_transid_exonobject[gene_name].setdefault(gene_id, dict())
        if current_ensembl_transcript not in gencode_genename_geneid_transid_exonobject[gene_name][gene_id]:
            gencode_genename_geneid_transid_exonobject[gene_name][gene_id].setdefault(current_ensembl_transcript, list())
            gencode_genename_geneid_transid_exonobject[gene_name][gene_id][current_ensembl_transcript].append(current_ensembl_exon)
        else:
            gencode_genename_geneid_transid_exonobject[gene_name][gene_id][current_ensembl_transcript].append(current_ensembl_exon)

    return (gencode_genename_geneid_transid_exonobject, current_ensembl_exon)


def parse_cds_gff3_line(
    feature_dict,
    current_ensembl_transcript,
    gencode_genename_geneid_transid_cdsobject,
    genome_v
):
    for inf in feature_dict["info"]:
        field = inf.split("=")[0].lower()
        result = inf.split("=")[1].upper()
        transcript_id = None

        if field == "gene_name":
            gene_name = result
        elif field == "id":
            result.replace("CDS:", "")
            transcript_id = result
            feature_dict["transcript_id"] = transcript_id

        elif field == "gene_id":
            gene_id = result

        feature_dict[field] = result
        feature_dict["genome_version"] = genome_v
    current_ensembl_cds = EnsemblCds(**{k: v for k, v in feature_dict.items() if k in EnsemblCds.__annotations__})

    if current_ensembl_cds is not None:
        current_ensembl_cds.get_coordinate_liftover()
        if gene_name not in gencode_genename_geneid_transid_cdsobject:
            gencode_genename_geneid_transid_cdsobject.setdefault(gene_name, dict())
        if gene_id not in gencode_genename_geneid_transid_cdsobject[gene_name]:
            gencode_genename_geneid_transid_cdsobject[gene_name].setdefault(gene_id, dict())
        if transcript_id not in gencode_genename_geneid_transid_cdsobject[gene_name][gene_id]:
            gencode_genename_geneid_transid_cdsobject[gene_name][gene_id].setdefault(current_ensembl_transcript, list())
            gencode_genename_geneid_transid_cdsobject[gene_name][gene_id][current_ensembl_transcript].append(current_ensembl_cds)
        else:
            gencode_genename_geneid_transid_cdsobject[gene_name][gene_id][current_ensembl_transcript].append(current_ensembl_cds)

    return (gencode_genename_geneid_transid_cdsobject, current_ensembl_cds)


def parse_gencode(
        gencode_gff: str,
        mane_clin_ens_id: set,
        mane_select_ens_id: set,
        lrg_ensembl: dict,
        genome_v: int
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
                    gencode_genename_geneobject,
                    genome_v
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
                    lrg_ensembl,
                    genome_v
                )
                # adding the transcript object to the gene obj
                current_ensembl_gene.add_transcript(current_ensembl_transcript)
            # parse exon and create EnsemblExon dataclass with its information
            elif feature_dict["feature"] == "exon":
                current_ensembl_transcript.increment_numb_exons()
                (
                    gencode_genename_geneid_transid_exonobject,
                    current_ensembl_exon
                ) = parse_exon_gff3_line(
                        feature_dict,
                        current_ensembl_transcript,
                        gencode_genename_geneid_transid_exonobject,
                        genome_v
                    )
                # adding exon obj to the transcript obj
                current_ensembl_transcript.add_exon(current_ensembl_exon)

            # parse cds and create EnsemblCds dataclass with its information
            elif feature_dict["feature"] == "cds":
                current_ensembl_transcript.increment_numb_cds()
                (
                    gencode_genename_geneid_transid_cdsobject,
                    current_ensembl_cds
                ) = parse_cds_gff3_line(
                    feature_dict,
                    current_ensembl_transcript,
                    gencode_genename_geneid_transid_cdsobject,
                    genome_v
                )
                # adding cds obj to the transcript obj
                current_ensembl_transcript.add_cds(current_ensembl_cds)
            
            # As in the gencode gff we have 3', start codon... features,
            # we don't set to None any object if a feature that we haven't
            #  parsed exists. (different to refseq)
    logging.info(f"{gencode_gff} has been parsed successfully")

    return (
        gencode_genename_geneobject,
        gencode_genename_geneid_transobject,
        gencode_genename_geneid_transid_exonobject,
        gencode_genename_geneid_transid_cdsobject
    )


if "__main__" == __name__:
    (
        mane_clin_ens_id,
        mane_select_ens_id,
        mane_clin_refseq_id,
        mane_select_refseq_id
    ) = get_mane(mane_gff)

    lrg_ensembl, lrg_refseq = get_lrg_trancripts(lrg_gff)

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
