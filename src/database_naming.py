import os
import sys
import gzip
from collections import defaultdict
from collections import defaultdict
from sqlalchemy import create_engine, exc
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy_utils import database_exists, create_database
from sqlalchemy import Column, Integer, String, Float, ForeignKey
import logging

# Configure the root logger
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Create a StreamHandler and set its output to sys.stdout
console_handler = logging.StreamHandler(sys.stdout)

# Optionally, you can set the logging level for the handler
console_handler.setLevel(logging.DEBUG)

# Add the handler to the root logger
logging.getLogger().addHandler(console_handler)

gene_table_cols = [
    "name",
    "id",
    "chr",
    "start",
    "end",
    "hgnc",
    "type"
]

gencode_gff = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/gencode/gencode.v43.annotation.gff3.gz"


def parse_gencode(gencode_gff):

    # Creates a dictionary with a default value of an empty dictionary
    gencode_gene_dict = defaultdict(dict)
    gencode_transcript_dict = defaultdict(dict)
    gencode_exon_dict = defaultdict(dict)
    gencode_cds_dict = defaultdict(dict)
    global_gene_name_dict = defaultdict(dict)

    with gzip.open(gencode_gff, 'rt') as fin:
        feature_fields = set()
        for line in fin:
            line = line.rstrip("\n")

            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            feature = tmp[2]
            chr = tmp[0].replace("chr", "")
            pos = tmp[3]
            end = tmp[4]
            info = tmp[8].split(";")

            if feature == "gene":
                gene_dict = dict()
                gene_dict["chromosome"] = chr
                gene_dict["start"] = pos
                gene_dict["end"] = end
                coordinates = f"{chr}:{pos}:{end}"

                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    feature_fields.add(field)

                    # print(result)
                    # if field == "gene_name":
                    #     if not result in gene_name:
                    #         gene_name.append(result)
                    #     else:
                    #         same_gene_name += 1

                    if field == "gene_name":
                        gene_name = result

                    elif field == "gene_id":
                        gene_id = result

                    # if field == "ID":
                        # there are genes with same identifier but with _PAR_Y
                        # at the final of the ID e.g. ENSG00000229232.6_PAR_Y
                        # (that we will differenciate between them as they are
                        # pseudoautosomal region on the Y chromosome
                        # if "PAR_Y" not in result.split(".")[1]:
                        #     # not taking into account gene version
                        #     result = result.split(".")[0]
                        # else:
                        #     result = str(result.split(".")[0]) + "_PAR_Y"
                    gene_dict[field] = result
                
                global_gene_name_dict["gencode"][gene_name]["gene"] = coordinates
                gencode_gene_dict[gene_id] = gene_dict

            elif feature == "transcript":
                transcript_dict = dict()
                transcript_dict["chromosome"] = chr
                transcript_dict["start"] = pos
                transcript_dict["end"] = end
                coordinates = f"{chr}:{pos}:{end}"

                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    feature_fields.add(field)

                    if field == "transcript_id":
                        transcirpt_id = result

                    elif field == "gene_name":
                        gene_name = result

                    transcript_dict[field] = result
                
                global_gene_name_dict["gencode"][gene_name]["transcript"] = coordinates
                gencode_transcript_dict[transcirpt_id] = transcript_dict

            elif feature == "exon":
                exon_dict = dict()
                exon_dict["chromosome"] = chr
                exon_dict["start"] = pos
                exon_dict["end"] = end
                coordinates = f"{chr}:{pos}:{end}"

                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    feature_fields.add(field)

                    if field == "exon_id":
                        exon_id = result

                    elif field == "gene_name":
                        gene_name = result

                    exon_dict[field] = result

                global_gene_name_dict["gencode"][gene_name]["exon"] = coordinates
                gencode_exon_dict[exon_id] = exon_dict

            elif feature == "CDS":
                cds_dict = dict()
                cds_dict["chromosome"] = chr
                cds_dict["start"] = pos
                cds_dict["end"] = end
                coordinates = f"{chr}:{pos}:{end}"

                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    feature_fields.add(field)

                    if field == "id":
                        exon_id = result

                    elif field == "gene_name":
                        gene_name = result

                    cds_dict[field] = result
                
                global_gene_name_dict["gencode"][gene_name]["cds"] = coordinates
                gencode_cds_dict[exon_id] = cds_dict
            

        
    return (gencode_gene_dict, gencode_transcript_dict, gencode_exon_dict, gencode_cds_dict)


refseq_gff = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/db_uri/refseq/GRCh38_latest_genomic.gff.gz"
def get_refseq_transcripts(refseq_gff: str, global_gene_name_dict):
    '''
        Function that reads a RefSeq gff3 file and returns two structured dictionaries
            with genes and transcripts
        :param str mane_gff: RefSeq gff3 file
        :return: two dicts, one with primary keys pointing to an enst_id, and a second one
            with gene id's as primary keys
        :rtype: Two defaultdict(dict)
    '''
    gene_numb = 0
    hgnc_numb = 0
    feature_set = set()
    refseq_transcripts = defaultdict(dict)
    refseq_genes_dict  = defaultdict(dict)
    refseq_cds_dict = defaultdict(dict)
    feature_dict = dict()
    with gzip.open(refseq_gff, 'rt') as fin:
        for line in fin:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            feature = tmp[2].lower()
            chr = tmp[0]
            pos = tmp[3]
            end = tmp[4]
            info= tmp[8].split(";")
            feature_dict["chromosome"] = chr
            feature_dict["start"] = pos
            feature_dict["end"] = end
            coordinates = f"{chr}:{pos}:{end}"
            if feature not in feature_set:
                feature_set.add(feature)

            if feature == "gene":
                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    if field == "dbxref":
                        dbxref_results = result.split(",")
                        for dbxref_result in dbxref_results:
                            if "hgnc" in dbxref_result:
                                hgnc_id = dbxref_result.split(":")[-1]
                                hgnc_id = f"HGNC:{hgnc_id}"
                                feature_dict["hgnc_id"] = hgnc_id
                                hgnc_numb += 1
                            elif "GeneID" in dbxref_result:
                                gene_id = dbxref_result.split(":")[1]
                                feature_dict["id"] = gene_id
                    elif field == "name":
                        gene_name = result
                    

                global_gene_name_dict[gene_name]["refseq"][gene_id] = coordinates

            if feature == "transcript" or feature == "primary_transcript" or feature == "lnc_rna" or feature == "mrna" or feature == "snrna":
                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    if field == "parent":
                        # removing gene- from the gene name
                        gene_name = result.split("-", 1)[-1]
                    if field == "id":
                        transcript_id = result
                    if field == "dbxref":
                        dbxref_results = result.split(",")
                        for dbxref_result in dbxref_results:
                            if "geneid" in dbxref_result:
                                gene_id = dbxref_result.split(":")[1]
                                feature_dict["id"] = gene_id
                            
                            # elif "hgnc" in dbxref_result:
                            #     hgnc_id = dbxref_result.split(":")[-1]
                            #     hgnc_id = f"HGNC:{hgnc_id}"
                            #     feature_dict["hgnc_id"] = hgnc_id
                            #     hgnc_numb += 1
                global_gene_name_dict[gene_name]["refseq"]["transcripts"][transcript_id]["coords"] = coordinates
                
            if feature == "exon":

                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    if field == "id":
                        exon_id = result
                    elif field == "gene":
                        gene_name = result
                    elif field == "parent":
                        transcript_id = result

                    # if field == "dbxref":
                    #     dbxref_results = result.split(",")
                    #     for dbxref_result in dbxref_results:
                    #         if "hgnc" in dbxref_result:
                    #             hgnc_id = dbxref_result.split(":")[-1]
                    #             hgnc_id = f"HGNC:{hgnc_id}"
                    #             feature_dict["hgnc_id"] = hgnc_id
                    #             hgnc_numb += 1
                    #         elif "geneid" in dbxref_result:
                    #             gene_id = dbxref_result.split(":")[1]
                    #             feature_dict["id"] = gene_id
                            # elif "mirbase" in dbxref_result:
                            #     mirbase_id = dbxref_result.split(":")[1]
                            #     feature_dict["mirbase_id"] = mirbase_id
                            # elif "mim" in dbxref_result:
                            #     mim_id = dbxref_result.split(":")[1]
                            #     feature_dict["mim_id"] = mim_id
                    elif field == "gene_biotype":
                        field = "gene_type"
                    elif field == "gene":
                        id = result
                    feature_dict[field] = result
                global_gene_name_dict[gene_name]["refseq"]["transcripts"][transcript_id]["exons"][exon_id] = coordinates

                # refseq_genes_dict[id] = feature_dict
            
            
            if feature == "cds":

                for inf in info:
                    field = inf.split("=")[0].lower()
                    result = inf.split("=")[1]
                    if field == "id":
                        cds_id = result
                    if field == "gene":
                        gene_name = result
                    if field == "parent":
                        transcript_id = result
                global_gene_name_dict[gene_name]["refseq"]["transcripts"][transcript_id]["cds"][cds_id] = coordinates
                
    print(global_gene_name_dict["MIR6859-1"].keys())
    print(global_gene_name_dict["MIR6859-1"]["refseq"].keys())
    print(global_gene_name_dict["MIR6859-1"]["refseq"]["transcripts"].keys())



    return refseq_genes_dict

# 
#   
# (
#     gencode_gene_dict,
#     gencode_transcript_dict,
#     gencode_exon_dict,
#     gencode_cds_dict
# ) = parse_gencode(gencode_gff)
global_gene_name_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))))
get_refseq_transcripts(refseq_gff, global_gene_name_dict)


def create_database():
    output_dir = "/home/ocanal/Desktop/gene_isoforms/build_gene_transcript_db/"
    mysql_psw = os.environ["MYSQL_PWD"]
    # Now define and create the genes database
    genes_db= output_dir + "/genes.db"
    engine  = create_engine(f'mysql+mysqlconnector://ocanal:{mysql_psw}@localhost:3306/genes_db')
    Session = sessionmaker(bind=engine)
    session = Session()
    Base = declarative_base()

    try:
        engine.execute("CREATE DATABASE genes_db")
        print("database genes_db created successfully")
    
    except exc.SQLAlchemyError as e:
        print("Error ocurred while creating the database:", str(e))

    class Genes(Base):
        __tablename__ = "Genes"
        gene_id = Column(String(50), primary_key=True)
        start = Column(Integer, nullable=False)
        end = Column(Integer, nullable=False)
        chromosome = Column(String(20), nullable=False)
        id = Column(String(50), nullable=False)
        tag = Column(String(100))
        gene_name = Column(String(50))
        gene_type = Column(String(50))
        level = Column(Integer)
        havana_gene = Column(String(50))
        hgnc_id = Column(String(50))

        transcripts = relationship("Transcripts", back_populates="gene")

        def __repr__(self):
            return f"Gene {self.gene_id}"

    class Transcripts(Base):
        __tablename__ = "Transcripts"
        transcript_id = Column(String(50), primary_key=True)
        start = Column(Integer, nullable=False)
        end = Column(Integer, nullable=False)
        chromosome = Column(String(10), nullable=False)
        transcript_type = Column(String(50))
        transcript_name = Column(String(50))
        id = Column(String(50), nullable=False)
        havana_transcript = Column(String(50))
        protein_id = Column(String(50))
        gene_id = Column(String(50), ForeignKey("Genes.gene_id"))

        # Relationships
        gene = relationship("Genes", back_populates="transcripts")
        exons = relationship("Exons", back_populates="transcript")

        def __repr__(self):
            return f"Transcript {self.transcript_id}"

    class Exons(Base):
        __tablename__ = "Exons"
        exon_id = Column(String(50), primary_key=True)
        id = Column(String(50), nullable=False)
        exon_number = Column(Integer)
        start = Column(Integer, nullable=False)
        end = Column(Integer, nullable=False)
        chromosome = Column(String(50), nullable=False)
        transcript_id = Column(String(50), ForeignKey("Transcripts.transcript_id"))

        #Relationships
        transcript = relationship("Transcripts", back_populates="exons")
        cds = relationship("Cds", back_populates="exon")

        def __repr__(self):
            return f"Exon {self.exon_id}"

    class Cds(Base):
        __tablename__ = "Cds"
        id = Column(String(50), primary_key=True)
        start = Column(Integer, nullable=False)
        end = Column(Integer, nullable=False)
        chromosome = Column(String(50), nullable=False)
        exon_id = Column(String(50), ForeignKey("Exons.exon_id"))

        # Relationships
        exon = relationship("Exons", back_populates="cds")

        def __repr__(self):
            return f"CDS: {self.id}"

    Base.metadata.create_all(engine)

    return Genes, Transcripts, Exons, Cds, session


# Genes, Transcripts, Exons, Cds, session = create_database()


def insert_db_gencode_gene_data(gencode_gene_dict, session, Genes):
    """
    From a dictionary with [gene_id]:[feature]:[values]
    inserts the features and values into the genes table

    Params:
        gencode_gene_dict[gene_id]:[feature]:[values]
        session: session of the database, returned by function create_database
    """
    valid_gene_keys = {
        'gene_id',
        'start',
        'end',
        'chromosome',
        'id',
        'tag',
        'gene_name',
        'gene_type',
        'level',
        "havana_gene",
        "hgnc_id"
    }

    for gene, gene_items in gencode_gene_dict.items():
        filtered_gene_data = {}
        
        # Filtering the dict to only contain the key-value pairs that are contained in the database
        filtered_gene_data = {key: value for key, value in gene_items.items() if key in valid_gene_keys}

        # Using ** it recieves the following data: gene_id="gene_id", start=100, end=200 ...
        new_gene = Genes(**filtered_gene_data)

        try:
            print(f"inserting gene {new_gene}")
            # Add new gene to the session
            session.add(new_gene)
            session.commit()
        except:
            msg = "Error in database writing"
            raise SQLAlchemyError(msg)
        else:
            msg = f"INFO: gene {new_gene} record successfully"
            logging.info(msg)


def insert_db_gencode_transcirpts_data(gencode_transcript_dict, session, Transcripts):
    """
    From a dictionary with [transcript_id]:[feature]:[values]
    inserts the features and values into the transcripts table

    Params:
        gencode_transcript_dict[transcript_id]:[feature]:[values]
        session: session of the database, returned by function create_database
        Transcripts: Transcripts database table
    """

    valid_transcript_keys = {
        'transcript_id',
        'start',
        'end',
        'chromosome',
        'transcript_type',
        'transcript_name',
        'id',
        'havana_transcript',
        'protein_id',
        "gene_id"
    }

    for transcirpt, transcript_items in gencode_transcript_dict.items():
        filtered_transcript_data = {}
        
        # Filtering the dict to only contain the key-value pairs that are contained in the database
        filtered_transcript_data = {key: value for key, value in transcript_items.items() if key in valid_transcript_keys}

        # Using ** it recieves the following data: gene_id="gene_id", start=100, end=200 ...
        new_transcript = Transcripts(**filtered_transcript_data)

        try:
            print(f"inserting transcript {new_transcript}")
            # Add new gene to the session
            session.add(new_transcript)
            session.commit()
        except:
            msg = "Error in database writing"
            raise SQLAlchemyError(msg)
        else:
            msg = f"INFO: gene {new_transcript} record successfully"
            logging.info(msg)

def insert_db_gencode_exons_data(gencode_exon_dict, session, Exons):
    """
    From a dictionary with [exon_id]:[feature]:[values]
    inserts the features and values into the Exons table

    Params:
        gencode_exon_dict[exon_id]:[feature]:[values]
        session: session of the database, returned by function create_database
        Exons: Exons database table
    """


    valid_exon_keys = {
        'exon_id',
        'start',
        'end',
        'chromosome',
        'exon_number',
        'id',
        'transcript_id'
    }

    for expon, exon_items in gencode_exon_dict.items():
        filtered_exon_data = {}
        
        # Filtering the dict to only contain the key-value pairs that are contained in the database
        filtered_exon_data = {key: value for key, value in exon_items.items() if key in valid_exon_keys}

        # Using ** it recieves the following data: gene_id="gene_id", start=100, end=200 ...
        new_exon = Exons(**filtered_exon_data)

        try:
            print(f"inserting exon {new_exon}")
            # Add new gene to the session
            session.add(new_exon)
            session.commit()
        except:
            msg = "Error in database writing"
            raise SQLAlchemyError(msg)
        else:
            msg = f"INFO: gene {new_exon} record successfully"
            logging.info(msg)


def insert_db_gencode_cds_data(gencode_cds_dict, session, Cds):
    """
    From a dictionary with [cds_id]:[feature]:[values]
    inserts the features and values into the Exons table

    Params:
        gencode_cds_dict[cds_id]:[feature]:[values]
        session: session of the database, returned by function create_database
        Cds: Cds database table
    """
    class Cds(Base):
        __tablename__ = "Cds"
        id = Column(String, primary_key=True)
        start = Column(Integer, nullable=False)
        end = Column(Integer, nullable=False)
        chromosome = Column(String, nullable=False)
        exon_id = Column(String, ForeignKey("Exons.exon_id"))

        # Relationships

    valid_cds_keys = {
        'id',
        'start',
        'end',
        'chromosome',
        'exon_id'
    }

    for cds, cds_items in gencode_cds_dict.items():
        filtered_cds_data = {}
        
        # Filtering the dict to only contain the key-value pairs that are contained in the database
        filtered_cds_data = {key: value for key, value in cds_items.items() if key in valid_cds_keys}

        # Using ** it recieves the following data: gene_id="gene_id", start=100, end=200 ...
        new_cds = Cds(**filtered_cds_data)

        try:
            print(f"inserting cds {new_cds}")
            # Add new gene to the session
            session.add(new_cds)
            session.commit()
        except:
            msg = "Error in database writing"
            raise SQLAlchemyError(msg)
        else:
            msg = f"INFO: gene {new_cds} record successfully"
            logging.info(msg)

# insert_db_gencode_gene_data(gencode_gene_dict, session, Genes)
# insert_db_gencode_transcirpts_data(gencode_transcript_dict, session, Transcripts)
# insert_db_gencode_exons_data(gencode_exon_dict, session, Exons)
# insert_db_gencode_cds_data(gencode_cds_dict, session, Cds)