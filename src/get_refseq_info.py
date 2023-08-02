from dataclasses import dataclass, field
import gzip

from src.global_variables import (
    convert_chromosomes,
    logging,
    check_genome_version,
    lo_38_to_37,
    lo_37_to_38
)
from src.get_mane import get_mane, mane_gff
from src.get_lrg_info import get_lrg_trancripts, lrg_gff

refseq_gff_grch38 = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/refseq/GRCh38_latest_genomic.gff.gz"
refseq_gff_grch37 = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/refseq/GRCh37_latest_genomic.gff.gz"
mock_refseq_gff_grch38 = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/refseq/mock.GRCh38_latest_genomic.gff.gz"
mock_refseq_gff_grch37 = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/refseq/mock.GRCh37_latest_genomic.gff.gz"

@dataclass
class RefseqGene:
    gene_id:                    str
    start:                      int
    end:                        int
    chromosome:                 int
    id:                         str
    name:                       str
    gene_type:                  str
    feature:                    str
    genome_version:             int
    last_assigned_primary_key:  int = 0
    gene_synonyms:              str = field(default=None)
    hgnc_id:                    int = field(default=None)
    num_transcripts:            int = field(default=0)
    gencode_same_gene_coords:   str = field(default=None)
    liftover_chromosome:        int = field(default=None)
    liftover_start:             int = field(default=None)
    liftover_end:               int = field(default=None)

    def __post_init__(self):
        RefseqGene.last_assigned_primary_key +=1
        self.primary_key = RefseqGene.last_assigned_primary_key

    def get_full_id(self):
        return (f"{self.id}_{self.genome_version}")
    
    def __hash__(self):
        return hash(self.gene_id)

    def __eq__(self, other):
        if isinstance(other, RefseqGene):
            return self.gene_id == other.gene_id
        return False

    def get_coords(self):
        return f"{self.chromosome}:{int(self.start)}-{int(self.end)}"

    def increment_num_transcripts(self):
        self.num_transcripts += 1

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
            print(self.chromosome, self.start, self.end)
            start = lo_37_to_38.convert_coordinate(
                chromosome_name,
                int(self.start)
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
# class RefseqGene_Grch37:
#     gene_id:                    str
#     start:                      int
#     end:                        int
#     chromosome:                 int
#     id:                         str
#     name:                       str
#     gene_type:                  str
#     feature:                    str
#     primary_key:                int = field(default=None)
#     gene_synonyms:              str = field(default=None)
#     hgnc_id:                    int = field(default=None)
#     num_transcripts:            int = field(default=0)
#     gencode_same_gene_coords:   str = field(default=None)
#     grch38_lo_chromosome:        int = field(default=None)
#     grch38_lo_start:             int = field(default=None)
#     grch38_lo_end:               int = field(default=None)
    

#     def __hash__(self):
#         return hash(self.gene_id)
    
#     def __eq__(self, other):
#         if isinstance(other, RefseqGene):
#             return self.gene_id == other.gene_id
#         return False

#     def get_coords(self):
#         return f"{self.chromosome}:{int(self.start)}-{int(self.end)}"

#     def increment_num_transcripts(self):
#         self.num_transcripts += 1

#     def get_coordinate_liftover(self):
#         grch37 = lo_37_to_38.convert_coordinate(
#             self.chromosome,
#             int(self.start),
#             int(self.end)
#         )
#         if grch37:
#             # liftover performed successfully
#             (
#                 self.grch38_lo_chromosome,
#                 self.grch38_lo_start,
#                 self.grch38_lo_end,
#                 _
#             ) = grch37[0]
#             logging.info(f"liftover performed successfully for {self.gene_id}")


@dataclass
class RefseqTranscript:
    id:                         str
    start:                      int
    end:                        int
    chromosome:                 str
    gene_id:                    int
    gene_name:                  str
    genome_version:             int
    parent:                     str
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
    liftover_chromosome:        int = field(default=None)
    liftover_start:             int = field(default=None)
    liftover_end:               int = field(default=None)

    def get_full_id(self):
        return (f"{self.id}_{self.genome_version}")
    
    def increment_numb_exons(self):
        self.numb_exons += 1

    def increment_numb_cds(self):
        self.numb_cds += 1

    def get_coords(self):
        return f"{self.chromosome}:{int(self.start)}-{int(self.end)}"

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        if isinstance(other, RefseqTranscript):
            return self.id == other.id
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
                    {self.liftover_chromosome}:{self.liftover_start}-{self.liftover_end}")
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
                    {self.liftover_chromosome}:{self.liftover_start}-{self.liftover_end}")
# @dataclass
# class RefseqTranscript_Grch37:
#     id:                         str
#     start:                      int
#     end:                        int
#     chromosome:                 str
#     gene_id:                    int
#     gene_name:                  str
#     hgnc_id:                    str = field(default=None)
#     lrg_id:                     str = field(default=None)
#     lrg_transcript:             str = field(default=None)
#     ccds:                       str = field(default=None)
#     numb_exons:                 int = field(default=0)
#     numb_cds:                   int = field(default=0)
#     mane_clin:                  bool = field(default=False)
#     mane_select:                bool = field(default=False)
#     gencode_same_trans_coords:  str = field(default=None)
#     gencode_same_cds_coords:    str = field(default=None)
#     grch38_lo_chromosome:       int = field(default=None)
#     grch38_lo_start:            int = field(default=None)
#     grch38_lo_end:              int = field(default=None)

#     def increment_numb_exons(self):
#         self.numb_exons += 1

#     def increment_numb_cds(self):
#         self.numb_cds += 1

#     def get_coords(self):
#         return f"{self.chromosome}:{int(self.start)}-{int(self.end)}"

#     def __hash__(self):
#         return hash(self.id)

#     def __eq__(self, other):
#         if isinstance(other, RefseqTranscript):
#             return self.id == other.id
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
class RefseqExons:
    exon_id:                str
    transcript_id:          str
    start:                  int
    end:                    int
    chromosome:             str
    exon_number:            int
    genome_version:         int
    liftover_chromosome:    int = field(default=None)
    liftover_start:         int = field(default=None)
    liftover_end:           int = field(default=None)

    def get_full_id(self):
        return (f"{self.exon_id}_{self.genome_version}")

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
# class RefseqExons_Grch37:
#     exon_id:                str
#     transcript_id:          str
#     start:                  int
#     end:                    int
#     chromosome:             str
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



@dataclass
class RefseqCds:
    id:                         str
    start:                      int
    end:                        int
    chromosome:                 str
    transcript_id:              str
    genome_version:             int
    last_assigned_primary_key:  int = 0
    liftover_chromosome:        int = field(default=None)
    liftover_start:             int = field(default=None)
    liftover_end:               int = field(default=None)

    def __post_init__(self):
        RefseqCds.last_assigned_primary_key += 1
        self.primary_key = RefseqCds.last_assigned_primary_key
    
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


def parse_gff3_common_fields(line, feature_dict):
    """
    Parsing the gff3 file and extracting the common fields
    Params:
    ------
        line: each line of the gff3 file
        feature_dict: dict where will be stored all the information
            of the corresponding line. In this case initially is empty
    Return:
    -------
        feature_dict: dict where information of the element is stored
    """
    line = line.rstrip("\n")
    if line.startswith("#"):
        return (None)
    tmp = line.split("\t")
    feature_dict["feature"] = tmp[2].lower()
    refseq_chr = tmp[0]
    refseq_chr = refseq_chr.split(".")[0]
    feature_dict["chromosome"] = convert_chromosomes(refseq_chr)
    feature_dict["start"] = tmp[3]
    feature_dict["end"] = tmp[4]
    feature_dict["info"] = tmp[8].split(";")

    return (feature_dict)


def parse_gene_lines(
    feature_dict,
    refseq_genename_geneobject,
    genome_v
):
    """
    Parsing gene lines and creating a RefseqGene object that will
    be added to the refseq_genename_geneobject dict
    
    Params:
    -------
        feature_dict: dict where all gene info is stored
        resfeq_genename_geneobject: dict where genes can be identified
            by its genenames
    
    Return:
    -------
        current_gene_object: RefseqGene instance that contains all the information
            of the gene being parsed
        refseq_genename_geneobject: dictionary with the current gene line added
            structure: [genename] = [list of genes objects]
    """
    gene_id = None
    for inf in feature_dict["info"]:
        field = inf.split("=")[0].lower()
        result = inf.split("=")[1].upper()

        if field == "dbxref":
            dbxref_results = result.split(",")
            for dbxref_result in dbxref_results:
                if "GENEID" in dbxref_result:
                    gene_id = dbxref_result.split(":")[1]
                    feature_dict["gene_id"] = gene_id
                elif "HGNC" in dbxref_result:
                    hgnc_id = dbxref_result.split(":")[-1]
                    hgnc_id = f"HGNC:{hgnc_id}"
                    feature_dict["hgnc_id"] = hgnc_id

        elif field == "gene_synonym":
            if "," in result:
                feature_dict["gene_synonyms"] = result.split(",")
            else:
                feature_dict["gene_synonyms"] = [result]
        elif field == "gene_biotype":
            field = "gene_type"
        feature_dict[field] = result

    if "name" not in feature_dict:
        feature_dict["name"] = feature_dict["feature"]

    if "gene_biotype" not in feature_dict:
        feature_dict["gene_type"] = feature_dict["feature"]

    # There is one gene that doesn't contain gene_id
    if "gene_id" not in feature_dict:
        return (None) # there was a continue before

    feature_dict["genome_version"] = genome_v
    # just include the key:value pairs in feature dict that are dfined in the dataclass RefseqGene
    current_gene_obj = RefseqGene(**{k: v for k, v in feature_dict.items() if k in RefseqGene.__annotations__})

    if current_gene_obj is not None:
        current_gene_obj.get_coordinate_liftover()

        if feature_dict["name"] not in refseq_genename_geneobject:
            gene_name = feature_dict["name"]
            refseq_genename_geneobject[gene_name] = list()
            refseq_genename_geneobject[gene_name].append(current_gene_obj)
        else:
            gene_name = feature_dict["name"]
            refseq_genename_geneobject[gene_name].append(current_gene_obj)

    return (current_gene_obj, refseq_genename_geneobject)


def parse_transcript_lines(
    feature_dict,
    mane_clin_refseq_id,
    mane_select_refseq_id,
    lrg_refseq,
    refseq_genename_geneid_transobject,
    current_gene_obj,
    genome_v
):
    """
    parsing transcript lines to extract information and create
    a RefseqTranscript instance
    Params:
    -------
        feature_dict: dict where the features of the transcripts will be stored
        mane_clin_refseq_id: list of refseq ids that are mane clin
        mane_select_refseq_id: list of refseq ids that are mane select
        refseq_genename_geneid_transobject: dictionary where transcript objects
            are stored in order to be identified by the same gene name
            structure: [genename][geneid]=[list of transcript objects]
    """
    for inf in feature_dict["info"]:
        field = inf.split("=")[0].lower()
        result = inf.split("=")[1].upper()
        if field == "parent":
            if "-" in result:
                # removing gene- from the gene name
                gene_name = result.split("-", 1)[-1]
                feature_dict["parent"] = gene_name
        elif field == "gbkey":
            feature_dict["transcript_type"] = result
        elif field == "gene":
            feature_dict["gene_name"] = result
        elif field == "id":
            # in V_gene_segment feature, it does not contain rna
            if "RNA-" not in result:
                result.replace("ID", "")
            result = result.replace("RNA-", "")
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
                if "GENEID" in dbxref_result:
                    gene_id = dbxref_result.split(":")[1]
                    feature_dict["gene_id"] = gene_id
                elif "HGNC" in dbxref_result:
                    hgnc_id = dbxref_result.split(":")[-1]
                    hgnc_id = f"HGNC:{hgnc_id}"
                    feature_dict["hgnc_id"] = hgnc_id
    if "name" not in feature_dict:
        feature_dict["name"] = feature_dict["feature"]
    if "gene_biotype" not in feature_dict:
        feature_dict["gene_biotype"] = feature_dict["feature"]

    # miRNA does not contain transcripts_ids so we will use it's id to identify them 
    if "transcript_id" not in feature_dict:
        feature_dict["transcript_id"] = id

    # There are 13 transcript that don't contain transcript_id, to solve it:
    # if "transcript_id" not in feature_dict:
    #     continue
    feature_dict["genome_version"] = genome_v
    # crete RefseqTranscript instance and adding it to the dict
    current_trans_obj = RefseqTranscript(**{k: v for k, v in feature_dict.items() if k in RefseqTranscript.__annotations__})

    if current_trans_obj is not None:
        current_trans_obj.get_coordinate_liftover()
        if gene_name not in refseq_genename_geneid_transobject:
            refseq_genename_geneid_transobject.setdefault(gene_name, dict())
        if current_gene_obj not in refseq_genename_geneid_transobject[gene_name]:
            refseq_genename_geneid_transobject[gene_name].setdefault(current_gene_obj, list())
            refseq_genename_geneid_transobject[gene_name][current_gene_obj].append(current_trans_obj)
        else:
            refseq_genename_geneid_transobject[gene_name][current_gene_obj].append(current_trans_obj)

    return (current_trans_obj, refseq_genename_geneid_transobject)


def parse_exon_lines(
    feature_dict,
    current_trans_obj,
    refseq_genename_geneid_transid_exonobject,
    genome_v
):
    """
    Parse exon lines of the refseq gff file inserting its
    information in the feature_dict, and adding the exon object to the
    refseq_genename_geneid_transid_exonobject dictionary.

    Params:
    -------
        feature_dict: dict where are stored information about the exon line
        current_trans_obj: transcript object that the exon that is being
            parsed belogs to.
        refseq_genename_geneid_transid_exonobject: dict where exon instance is
        stored. Structure: [genename][geneid][transid]=[list of exonobjects]

    Returns:
    --------
        refseq_genename_geneid_transid_exonobject: dict where the current exon
            is stored with the structure given in the previous line

    """
    current_exon_obj = None

    # feature_dict["info"] contains information as a list containing
    # field=result in each item of the list
    # [ID=gene-BPY2C, Dbxref=GeneID:442868,HGNC:HGNC:18225, Name=BPY2C, ...]
    for inf in feature_dict["info"]:
        field = inf.split("=")[0].lower()
        result = inf.split("=")[1].upper()

        if field == "id":
            if "EXON-" in result:
                exon_id = result.replace("EXON-", "")
            elif "ID-" in result:
                exon_id = result.replace("ID-", "")
            feature_dict["exon_id"] = exon_id
            if "-" in exon_id:
                last_hyphen_index = exon_id.rindex("-")
                feature_dict["exon_number"] = exon_id[last_hyphen_index + 1:]
            else:
                feature_dict["exon_number"] = None

        elif field == "gene":
            gene_name = result
        elif field == "transcript_id":
            transcript_id = result
            feature_dict["transcript_id"] = transcript_id
        elif field == "parent":
            # in V_gene_segment feature, it dgencode_gene_synonyms_associationin rna
            if "RNA-" not in result:
                parent = result.replace("ID-", "")
            parent = result.replace("RNA-", "")

        # as seen in the example, dbxref contain multiple fields seperated by ,
        elif field == "dbxref":
            dbxref_results = result.split(",")
            for dbxref_result in dbxref_results:
                if "GENEID" in dbxref_result:
                    gene_id = dbxref_result.split(":")[1]
                    feature_dict["id"] = gene_id

        elif field == "gene_biotype":
            field = "gene_type"

        # adding all the fields to the dictionary
        feature_dict[field] = result
    # if no transcript_id field, we use the parent field as transcript_id
    if "transcript_id" not in feature_dict:
        # In some exons the parent of the exon is a gene and we don't want
        # this reference. No transcript is associated. e.g.
        if "GENE-" not in parent:
            feature_dict["transcript_id"] = parent
        else:
            logging.info(
                f"Exon with parent {parent} does not contain a"
                " transcript id associated and the parent seems to be a gene")
            return (refseq_genename_geneid_transid_exonobject)
    feature_dict["genome_version"] = genome_v

    # Adding only feature dict key:values which keys are defined
    # in RefseqExons dataclass
    current_exon_obj = RefseqExons(**{k: v for k, v in feature_dict.items() if k in RefseqExons.__annotations__})

    # Creating refseq_genename_geneid_transid_exonobject dict
    if current_exon_obj is not None:
        current_exon_obj.get_coordinate_liftover()
        if gene_name not in refseq_genename_geneid_transid_exonobject:
            refseq_genename_geneid_transid_exonobject.setdefault(gene_name, dict())
        if gene_id not in refseq_genename_geneid_transid_exonobject[gene_name]:
            refseq_genename_geneid_transid_exonobject[gene_name].setdefault(gene_id, dict())
        if current_trans_obj not in refseq_genename_geneid_transid_exonobject[gene_name][gene_id]:
            refseq_genename_geneid_transid_exonobject[gene_name][gene_id].setdefault(current_trans_obj, list())
            refseq_genename_geneid_transid_exonobject[gene_name][gene_id][current_trans_obj].append(current_exon_obj)
        else:
            refseq_genename_geneid_transid_exonobject[gene_name][gene_id][current_trans_obj].append(current_exon_obj)

    return (
        refseq_genename_geneid_transid_exonobject
    )


def parse_cds_lines(
    feature_dict,
    current_trans_obj,
    refseq_genename_geneid_transid_cdsobject,
    genome_v
):
    for inf in feature_dict["info"]:
        field = inf.split("=")[0].lower()
        result = inf.split("=")[1].upper()
        if field == "id":
            feature_dict["id"] = result
        elif field == "gene":
            gene_name = result
        elif field == "parent":
            feature_dict["transcript_id"] = result.replace("RNA-", "")
        elif field == "dbxref":
            dbxref_results = result.split(",")
            for dbxref_result in dbxref_results:
                if "GENEID" in dbxref_result:
                    gene_id = dbxref_result.split(":")[1]
                    feature_dict["id"] = gene_id
    feature_dict[field] = result

    feature_dict["genome_version"] = genome_v
    current_cds_obj = RefseqCds(**{k: v for k, v in feature_dict.items() if k in RefseqCds.__annotations__})

    if current_cds_obj is not None:
        current_cds_obj.get_coordinate_liftover()
        if gene_name not in refseq_genename_geneid_transid_cdsobject:
            refseq_genename_geneid_transid_cdsobject.setdefault(gene_name, dict())
        if gene_id not in refseq_genename_geneid_transid_cdsobject[gene_name]:
            refseq_genename_geneid_transid_cdsobject[gene_name].setdefault(gene_id, dict())
        if current_trans_obj not in refseq_genename_geneid_transid_cdsobject[gene_name][gene_id]:
            refseq_genename_geneid_transid_cdsobject[gene_name][gene_id].setdefault(current_trans_obj, list())
            refseq_genename_geneid_transid_cdsobject[gene_name][gene_id][current_trans_obj].append(current_cds_obj)
        else:
            refseq_genename_geneid_transid_cdsobject[gene_name][gene_id][current_trans_obj].append(current_cds_obj)

    return (
        feature_dict,
        current_cds_obj,
        refseq_genename_geneid_transid_cdsobject
    )


def get_gene_synonyms(feature_dict, gene_name_to_synonym_dict):
    
    gene_name = feature_dict["name"]
    # Adding gene synonyms to the dict
    if "gene_synonyms" in feature_dict:
        gene_synonyms = set(feature_dict["gene_synonyms"])
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

    return (gene_name_to_synonym_dict)


def parse_refseq_gff(
    refseq_gff: str,
    mane_clin_refseq_id: list,
    mane_select_refseq_id: list,
    lrg_refseq: dict,
    genome_v: int,
    gene_name_to_synonym_dict: dict
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
    # checking that genome version is 37,19 or 38
    check_genome_version(genome_v)
    refseq_genename_geneobject = dict()
    refseq_genename_geneid_transobject = dict()
    refseq_genename_geneid_transid_exonobject = dict()
    refseq_genename_geneid_transid_cdsobject = dict()
    gene_name_to_synonym_dict = dict()
    current_gene_obj = None
    current_trans_obj = None

    genes_types = {
        'gene',
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
            feature_dict = dict()
            feature_dict = parse_gff3_common_fields(line, feature_dict)
            # skipping lines that start with #, where the function
            # will return None
            if feature_dict is None:
                continue

            # parsing gene info
            if feature_dict["feature"] in genes_types:
                # if line feature = gene, parsing it and adding its object to
                # refseq_genename_geneobject dictionary
                result = parse_gene_lines(
                    feature_dict,
                    refseq_genename_geneobject,
                    genome_v
                )

                if result is None:
                    logging.info(
                        f"{feature_dict['name']} does not contain an id and will be discarted"
                        )
                    current_gene_obj = None
                else:
                    (
                        current_gene_obj,
                        refseq_genename_geneobject
                    ) = result

                # creating gene_synonyms dict based on synonyms given in the
                # refseq gff3 file
                gene_name_to_synonym_dict = get_gene_synonyms(
                    feature_dict,
                    gene_name_to_synonym_dict
                )

            elif feature_dict["feature"] in transcripts_types:
                current_trans_obj = None
                # set at the beggining because we have to increment the number of exons and cds 
                # for each transcript and the object 
                if current_gene_obj is not None:
                    current_gene_obj.increment_num_transcripts()
                    (
                        current_trans_obj,
                        refseq_genename_geneid_transobject
                    ) = parse_transcript_lines(
                        feature_dict,
                        mane_clin_refseq_id,
                        mane_select_refseq_id,
                        lrg_refseq,
                        refseq_genename_geneid_transobject,
                        current_gene_obj,
                        genome_v
                    )
            elif feature_dict["feature"] == "exon":
                if current_gene_obj is None:
                    continue
                if current_trans_obj is not None:
                    current_trans_obj.increment_numb_exons()
                    (
                        refseq_genename_geneid_transid_exonobject
                    ) = parse_exon_lines(
                        feature_dict,
                        current_trans_obj,
                        refseq_genename_geneid_transid_exonobject,
                        genome_v
                    )
            elif feature_dict["feature"] == "cds":
                current_cds_obj = None
                if current_gene_obj is None:
                    continue
                if current_trans_obj is not None:
                    current_trans_obj.increment_numb_cds()
                    (
                        feature_dict,
                        current_cds_obj,
                        refseq_genename_geneid_transid_cdsobject
                    ) = parse_cds_lines(
                        feature_dict,
                        current_trans_obj,
                        refseq_genename_geneid_transid_cdsobject,
                        genome_v
                    )
            # if the element type is not gene_types, transcript_types, exon or cds
            # all current objects will be set to None. This avoids counting exons or cds
            # of a transcripts that does not correspond to
            else:
                current_gene_obj = None
                current_trans_obj = None

    return (
        refseq_genename_geneobject,
        refseq_genename_geneid_transobject,
        refseq_genename_geneid_transid_exonobject,
        refseq_genename_geneid_transid_cdsobject,
        gene_name_to_synonym_dict
    )


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
