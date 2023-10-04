import os

from sqlalchemy import (
    create_engine,
    exc,
    Column,
    Integer,
    String,
    ForeignKey,
    Boolean,
    Table,
    Index
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship, foreign
from sqlalchemy_utils import create_database

from src.global_variables import logging

genes_db_name = "Ensembl_Refseq_db"


def create_mysql_database():
    mysql_psw = os.environ["MYSQL_PWD"]
    # Now define and create the genes database

    logging.info(f"Creating database {genes_db_name}")

    engine  = create_engine(f'mysql+mysqlconnector://ocanal:{mysql_psw}@localhost:3306')
    Session = sessionmaker(bind=engine)
    session = Session()
    Base = declarative_base()

    try:
        engine.execute(f"CREATE DATABASE IF NOT EXISTS {genes_db_name}")
        engine.execute(f"USE {genes_db_name}")
        logging.info(f"database {genes_db_name} created successfully")

    except exc.SQLAlchemyError as e:
        logging.info("Error ocurred while creating the database:", str(e))

    gencode_gene_synonyms_association = Table(
        "gencode_gene_synonyms_association",
        Base.metadata,
        Column(
            "gencode_gene_id_genome_version",
            String(50),
            ForeignKey("gencode_genes.gene_id_genome_version")
        ),
        Column(
            "synonym",
            String(50),
            ForeignKey("gene_synonyms.synonym"),
        )
    )

    refseq_gene_synonyms_association = Table(
        "refseq_gene_synonyms_association",
        Base.metadata,
        Column(
            "refseq_gene_primary_key",
            Integer,
            ForeignKey("refseq_genes.primary_key"),
        ),
        Column(
            "synonym",
            String(50),
            ForeignKey("gene_synonyms.synonym"),
        )
    )

    class File_version(Base):
        __tablename__ = "file_versions"
        database = Column(String(50), primary_key=True)
        version = Column(String(50))
        file_path = Column(String(50))

    class Gene_synonyms(Base):
        __tablename__ = "gene_synonyms"
        synonym = Column(String(50), primary_key=True)

        # genes = relationship("Genes", secondary=gencode_gene_synonyms_association, backref="synonims")

    class Gencode_genes(Base):
        __tablename__ = "gencode_genes"
        gene_id_genome_version = Column(String(50), primary_key=True)
        id = Column(String(50), index=True)
        genome_version = Column(Integer)
        start_grch38 = Column(Integer)
        end_grch38 = Column(Integer)
        chromosome_grch38 = Column(String(50))
        chromosome_grch37 = Column(String(50))
        start_grch37 = Column(Integer)
        end_grch37 = Column(Integer)
        gene_id = Column(String(50), nullable=False)
        tag = Column(String(100))
        num_gencode_transcripts = Column(Integer)
        gene_name = Column(String(50), index=True)
        gene_type = Column(String(50))
        level = Column(Integer)
        havana_gene = Column(String(50))
        hgnc_id = Column(String(50))
        
        # gene_grch37 = relationship(
        #     "Genes_Grch37",
        #     uselist=False,
        #     back_populates="gene_grch38"
        # )

        # refseq_same_gene_coords = relationship(
        #     "RefseqGenes",
        #     secondary=Refseq_Gencode_genes_association,
        #     backref="gencode_same_gene_coords",
        #     overlaps="refseq_genes,refseq_same_gene_coords",
        #     foreign_keys=[
        #         Refseq_Gencode_genes_association.c.gencode_gene_id,
        #         Refseq_Gencode_genes_association.c.gencode_genome_version,
        #     ],
        # )

        synonyms = relationship(
            "Gene_synonyms",
            secondary=gencode_gene_synonyms_association,
            backref="gencode_associated_genes"
        )
        transcripts = relationship(
            "Gencode_transcripts",
            back_populates="gene"
        )

        def __repr__(self):
            return f"Gene {self.gene_id_genome_version}"


    class Refseq_genes(Base):
        __tablename__ = "refseq_genes"
        primary_key = Column(Integer, primary_key=True)
        genome_version = Column(Integer, nullable=False)
        gene_id = Column(String(50), nullable=False, index=True)
        gene_name = Column(String(50), nullable=False, index=True)
        chromosome_grch38 = Column(String(50))
        start_grch38 = Column(Integer)
        end_grch38 = Column(Integer)
        chromosome_grch37 = Column(String(50))
        start_grch37 = Column(Integer)
        end_grch37 = Column(Integer)
        gene_type = Column(String(50))
        feature = Column(String(50))
        num_transcripts = Column(Integer)
        name = Column(String(50))
        hgnc_id = Column(String(50))
        gene_biotype = Column(String(50))


        # Relationships
        # one to many with transcripts:
        transcripts = relationship(
            "Refseq_transcripts",
            back_populates="gene"
        )
        # Many to many relationships with synonyms and gencode genes
        synonyms = relationship(
            "Gene_synonyms",
            secondary=refseq_gene_synonyms_association,
            backref="genes_synonyms"
        )

        def __repr__(self):
            return (
                f"RefseqGene primary key :{self.primary_key} , gene_id:\
                {self.gene_id}")
        # gencode_same_gene_coords = relationship(
        #     "Genes",
        #     secondary=Refseq_Gencode_genes_association,
        #     backref="gencode_genes"
        # )

    # Association table between refseq and gencode transcripts
    association_table_trans_coords = Table(
        'transcripts_refseq_gencode_same_coords',
        Base.metadata,
        Column(
            'gencode_transcript_id',
            String(50),
            ForeignKey('gencode_transcripts.transcript_id_genome_version'),
        ),
        Column(
            'refseq_transcript_id',
            String(50),
            ForeignKey('refseq_transcripts.id_genome_version'),
        ),
    )
    association_table_exon_coords = Table(
        'transcripts_refseq_gencode_same_exon_coords',
        Base.metadata,
        Column(
            'gencode_transcript_id',
            String(50),
            ForeignKey('gencode_transcripts.transcript_id_genome_version'),
        ),
        Column(
            'refseq_transcript_id',
            String(50),
            ForeignKey('refseq_transcripts.id_genome_version'),
        ),

    )
    association_table_cds_coords = Table(
        'transcripts_refseq_gencode_same_cds_coords',
        Base.metadata,
        Column(
            'gencode_transcript_id',
            String(50),
            ForeignKey('gencode_transcripts.transcript_id_genome_version'),
        ),
        Column(
            'refseq_transcript_id',
            String(50),
            ForeignKey('refseq_transcripts.id_genome_version'),
        ),
    )

    class Refseq_transcripts(Base):
        __tablename__ = "refseq_transcripts"
        id_genome_version = Column(String(50), primary_key=True)
        id = Column(String(50), index=True)
        genome_version = Column(Integer)
        start_grch38 = Column(Integer)
        end_grch38 = Column(Integer)
        chromosome_grch38 = Column(String(50))
        chromosome_grch37 = Column(String(50))
        start_grch37 = Column(Integer)
        end_grch37 = Column(Integer)
        gene_id = Column(Integer, index=True)
        gene_name = Column(String(50), index=True)
        hgnc_id = Column(String(50))
        lrg_id = Column(String(50))
        lrg_transcript = Column(String(50))
        ccds = Column(String(50))
        numb_exons = Column(Integer)
        numb_cds = Column(Integer)
        mane_clin = Column(Boolean)
        mane_select = Column(Boolean)
        gencode_mane_select_id = Column(String(50))

        # Foreign key to establish relationship with RefseqGene
        gene_pk = Column(Integer, ForeignKey("refseq_genes.primary_key"))

        # refseq_transcript_grch37 = relationship(
        #     "RefseqTranscripts_Grch37",
        #     uselist=False,
        #     back_populates="refseq_transcript_grch38"
        # )
        # One to many relationships
        exons = relationship(
            "Refseq_exons", 
            back_populates="transcript"  # name of the exons column
        )
        cds = relationship(
            "Refseq_cds",
            back_populates="transcript"  # name of cds column
        )
        gene = relationship(
            "Refseq_genes",
            back_populates="transcripts"
        )

        # Many to many relationships
        gencode_same_trans_coords = relationship(
            "Gencode_transcripts",
            secondary=association_table_trans_coords,
            backref="gencode_transcripts",
        )
        gencode_same_exons_coords = relationship(
            "Gencode_transcripts",
            secondary=association_table_exon_coords,
            backref="gencode_trans_same_exons"
        )
        gencode_same_cds_coords = relationship(
            "Gencode_transcripts",
            secondary=association_table_cds_coords,
            backref="gencode_trans_same_cds"
        )
        
        def __repr__(self):
            return f"Refseq Transcript{self.id_genome_version}"

    class Refseq_exons(Base):
        __tablename__ = "refseq_exons"
        exon_id_genome_version = Column(String(50), primary_key=True)
        exon_id = Column(String(50), index=True)
        genome_version = Column(Integer)
        start_grch38 = Column(Integer)
        end_grch38 = Column(Integer)
        chromosome_grch38 = Column(String(50))
        chromosome_grch37 = Column(String(50))
        start_grch37 = Column(Integer)
        end_grch37 = Column(Integer)
        exon_number = Column(Integer)

        # Foreign key to relate with Transcript
        transcript_id = Column(
            String(50),
            ForeignKey("refseq_transcripts.id"),
            nullable=False,
            foreign_keys=exon_id
        )

        # refseq_exon_grch37 = relationship(
        #     "RefseqExons_Grch37",
        #     back_populates="refseq_exon_grch38",
        #     uselist=False
        # )
        # Relationships
        transcript = relationship(
            "Refseq_transcripts",
            back_populates="exons"
        )

        def __repr__(self):
            return f"Refseq Exon {self.exon_id_genome_version}"

    class Refseq_cds(Base):
        __tablename__ = "refseq_cds"
        primary_key = Column(Integer, primary_key=True)
        genome_version = Column(Integer)
        id = Column(String(50))
        start_grch38 = Column(Integer)
        end_grch38 = Column(Integer)
        chromosome_grch38 = Column(String(50))
        chromosome_grch37 = Column(String(50))
        start_grch37 = Column(Integer)
        end_grch37 = Column(Integer)

        # Foreign key to relate with Transcript
        transcript_id = Column(
            String(50),
            ForeignKey("refseq_transcripts.id")
        )

        # Relationships
        transcript = relationship(
            "Refseq_transcripts",
            back_populates="cds"
        )

        def __repr__(self):
            return f"Refseq Cds {self.primary_key}"

    class Gencode_transcripts(Base):
        __tablename__ = "gencode_transcripts"
        transcript_id_genome_version = Column(String(50), primary_key=True)
        transcript_id = Column(String(50), index=True)
        genome_version = Column(Integer)
        mane_select = Column(Boolean, nullable=True, default=None)
        refseq_mane_select_id = Column(String(50))
        mane_clin = Column(Boolean, nullable=True, default=None)
        start_grch38 = Column(Integer)
        end_grch38 = Column(Integer)
        chromosome_grch38 = Column(String(50))
        chromosome_grch37 = Column(String(50))
        start_grch37 = Column(Integer)
        end_grch37 = Column(Integer)
        transcript_type = Column(String(50))
        transcript_name = Column(String(50))
        id = Column(String(50), nullable=False)
        havana_transcript = Column(String(50))
        protein_id = Column(String(50))
        gene_id = Column(String(50), ForeignKey("gencode_genes.id"))
        lrg_transcript = Column(String(50))
        lrg_id = Column(String(50))
        ccds = Column(String(50))
        numb_exons = Column(Integer)
        numb_cds = Column(Integer)


        # transcript_grch37 = relationship(
        #     "Transcripts_Grch37",
        #     uselist=False,
        #     back_populates="transcript_grch38"
        # )
        # One to many relationships
        gene = relationship(
            "Gencode_genes",
            back_populates="transcripts"  # column name of Genes relationship
        )
        exons = relationship(
            "Gencode_exons",
            back_populates="transcript"  # column name of exons relationship
            )
        cds = relationship(
            "Gencode_cds",
            back_populates="transcript"  # column name  of cds relationship
        )
        # Many to many relationships
        refseq_same_trans_coords = relationship(
            "Refseq_transcripts",
            secondary=association_table_trans_coords,
            backref="refseq_transcripts_same_trans_coords",
        )
        refseq_same_exons_coords = relationship(
            "Refseq_transcripts",
            secondary=association_table_exon_coords,
            backref="refseq_trans_same_exons_coords"
        )
        refseq_same_cds_coords = relationship(
            "Refseq_transcripts",
            secondary=association_table_cds_coords,
            backref="refseq_trans_same_cds_coords",
            overlaps="gencode_same_cds_coords,gencode_trans_same_cds"
        )

        def __repr__(self):
            return f"Transcript {self.transcript_id_genome_version}"

    class Gencode_exons(Base):
        __tablename__ = "gencode_exons"
        exon_id_genome_version = Column(String(50), primary_key=True)
        exon_id = Column(String(50), index=True)
        genome_version = Column(Integer)
        id = Column(String(50), nullable=False)
        exon_number = Column(Integer)
        start_grch38 = Column(Integer)
        end_grch38 = Column(Integer)
        chromosome_grch38 = Column(String(50))
        chromosome_grch37 = Column(String(50))
        start_grch37 = Column(Integer)
        end_grch37 = Column(Integer)
        transcript_id = Column(
            String(50),
            ForeignKey("gencode_transcripts.transcript_id_genome_version"),
            nullable=True
        )

        # Relationships
        transcript = relationship(
            "Gencode_transcripts",
            back_populates="exons"
        )

        # exon_grch37 = relationship(
        #     "Exons_Grch37",
        #     uselist=False,
        #     back_populates="exon_grch38",
        # )

        def __repr__(self):
            return f"Exon {self.exon_id_genome_version}"

    class Gencode_cds(Base):
        __tablename__ = "gencode_cds"
        primary_key = Column(Integer, primary_key=True)
        genome_version = Column(Integer)
        id = Column(String(50), index=True)
        start_grch38 = Column(Integer)
        end_grch38 = Column(Integer)
        chromosome_grch38 = Column(String(50))
        chromosome_grch37 = Column(String(50))
        start_grch37 = Column(Integer)
        end_grch37 = Column(Integer)
        transcript_id = Column(
            String(50),
            ForeignKey("gencode_transcripts.transcript_id_genome_version")
        )
        exon_id = Column(String(50))
        exon_number = Column(Integer)

        # Relationships
        transcript = relationship(
            "Gencode_transcripts",
            back_populates="cds"
        )

        def __repr__(self):
            return f"CDS: {self.primary_key}"

    # gencode_37_gene_synonyms_association = Table(
    #     "gencode_grch37_gene_synonyms_association",
    #     Base.metadata,
    #     Column(
    #         "gencode_Grch37_gene_id",
    #         String(50),
    #         ForeignKey("Genes_Grch37.id"),
    #     ),
    #     Column(
    #         "synonym",
    #         String(50),
    #         ForeignKey("Gene_Synonyms.synonym"),
    #     )
    # )
    # refseq_37_gene_synonyms_association = Table(
    #     "refseq_grch37_gene_synonyms_association",
    #     Base.metadata,
    #     Column(
    #         "refseq_gene_primary_key",
    #         Integer,
    #         ForeignKey("RefseqGenes_Grch37.primary_key"),
    #     ),
    #     Column(
    #         "synonym",
    #         String(50),
    #         ForeignKey("Gene_Synonyms.synonym"),
    #     )
    # )
    # refseq_Gencode_37_genes_association = Table(
    #     "Refseq_Gencode_Grch37_genes_association",
    #     Base.metadata,
    #     Column(
    #         "gencode_Grch37_gene_id",
    #         String(50),
    #         ForeignKey("Genes_Grch37.id")
    #         ),
    #     Column(
    #         "refseq_gene_Grch37_primary_key",
    #         Integer,
    #         ForeignKey("RefseqGenes_Grch37.primary_key")
    #     )
    # )

    # class Genes_Grch37(Base):
    #     __tablename__ = "Genes_Grch37"
    #     id = Column(String(50), primary_key=True)
    #     start = Column(Integer, nullable=False)
    #     end = Column(Integer, nullable=False)
    #     chromosome = Column(String(20), nullable=False)
    #     gene_id = Column(String(50), nullable=False)
    #     tag = Column(String(100))
    #     num_gencode_transcripts = Column(Integer)
    #     num_refseq_transcripts = Column(Integer)
    #     gene_name = Column(String(50))
    #     gene_type = Column(String(50))
    #     level = Column(Integer)
    #     havana_gene = Column(String(50))
    #     hgnc_id = Column(String(50))
    #     gene_grch38_id = Column(String(50), ForeignKey("Genes.id"))
    #     gene_grch38 = relationship(
    #         "Genes",
    #         back_populates="gene_grch37",
    #     )

    #     refseq_same_gene_coords = relationship(
    #         "RefseqGenes_Grch37",
    #         secondary=refseq_Gencode_37_genes_association,
    #         backref="refseq_genes",
    #         overlaps="refseq_genes,refseq_same_gene_coords",
    #         secondaryjoin=(
    #             "and_("
    #             "Genes_Grch37.id == Refseq_Gencode_genes_association.c.gencode_Grch37_gene_id, "
    #             "RefseqGenes_Grch37.primary_key == Refseq_Gencode_genes_association.c.refseq_gene_Grch37_primary_key)"
    #         )
    #     )
    #     synonyms = relationship("Gene_Synonyms", secondary=gencode_37_gene_synonyms_association, backref="Genes_Grch37")
    #     transcripts = relationship("Transcripts", back_populates="gene")

    #     def __repr__(self):
    #         return f"Gene {self.gene_id}"
        
    # class RefseqGenes_Grch37(Base):
    #     __tablename__ = "RefseqGenes_Grch37"
    #     primary_key = Column(Integer, primary_key=True, autoincrement=True)
    #     gene_id = Column(String(50))
    #     gene_name = Column(String(50), nullable=False)
    #     start = Column(Integer, nullable=False)
    #     end = Column(Integer, nullable=False)
    #     chromosome = Column(String(50), nullable=False)
    #     gene_type = Column(String(50))
    #     feature = Column(String(50))
    #     num_transcripts = Column(Integer)
    #     name = Column(String(50))
    #     hgnc_id = Column(String(50))
    #     gene_biotype = Column(String(50))


    #     # Relationships

    #     # one to many with transcripts:
    #     transcripts = relationship(
    #         "RefseqTranscripts_Grch37",
    #         back_populates="gene"
    #     )
    #     # Many to many relationships with synonyms and gencode genes
    #     synonyms = relationship(
    #         "Gene_Synonyms",
    #         secondary=refseq_37_gene_synonyms_association,
    #         backref="genes_synonyms"
    #     )
    #     gencode_same_gene_coords = relationship(
    #         "Genes_Grch37",
    #         secondary=refseq_Gencode_37_genes_association,
    #         backref="gencode_genes",
    #         secondaryjoin=(
    #             "and_("
    #             "RefseqGenes_Grch37.gencode_Grch37_gene_id == "
    #             "Refseq_Gencode_Grch37_genes_association.c.gencode_Grch37_gene_id, "
    #             "RefseqGenes_Grch37.primary_key == "
    #             "Refseq_Gencode_Grch37_genes_association.c.refseq_gene_Grch37_primary_key)"
    #         )
    #     )

    # # Association table between refseq and gencode transcripts
    # association_table_trans_coords_Grch37 = Table(
    #     'transcripts_refseq_gencode_Grch37_same_coords',
    #     Base.metadata,
    #     Column(
    #         'transcript_id',
    #         String(50),
    #         ForeignKey('Transcripts_Grch37.transcript_id')
    #     ),
    #     Column(
    #         'id',
    #         String(50),
    #         ForeignKey('RefseqTranscripts_Grch37.id')
    #     )
    # )
    # association_table_exon_coords_Grch37 = Table(
    #     'transcripts_refseq_gencode_Grch37_same_exon_coords',
    #     Base.metadata,
    #     Column(
    #         'transcript_id',
    #         String(50),
    #         ForeignKey('Transcripts_Grch37.transcript_id')
    #     ),
    #     Column(
    #         'id',
    #         String(50),
    #         ForeignKey('RefseqTranscripts_Grch37.id')
    #     )
    # )
    # association_table_cds_coords_Grch37 = Table(
    #     'transcripts_refseq_gencode_Grch37_same_cds_coords',
    #     Base.metadata,
    #     Column(
    #         'transcript_id',
    #         String(50),
    #         ForeignKey('Transcripts_Grch37.transcript_id')
    #     ),
    #     Column(
    #         'id',
    #         String(50),
    #         ForeignKey('RefseqTranscripts_Grch37.id')
    #     )
    # )

    # class RefseqTranscripts_Grch37(Base):
    #     __tablename__ = "RefseqTranscripts_Grch37"
    #     id = Column(String(50), primary_key=True)
    #     start = Column(Integer, nullable=False)
    #     end = Column(Integer, nullable=False)
    #     chromosome = Column(String(50), nullable=False)
    #     gene_id = Column(Integer)
    #     gene_name = Column(String(50))
    #     hgnc_id = Column(String(50))
    #     lrg_id = Column(String(50))
    #     lrg_transcript = Column(String(50))
    #     ccds = Column(String(50))
    #     numb_exons = Column(Integer)
    #     numb_cds = Column(Integer)
    #     mane_clin = Column(Boolean)
    #     mane_select = Column(Boolean)
    #     # Foreign key to establish relationship with RefseqGene
    #     gene_pk = Column(Integer, ForeignKey("RefseqGenes_Grch37.primary_key"))
    #     refseq_transcript_grch38_id = Column(
    #         String(50),
    #         ForeignKey("RefseqTranscripts.id")
    #     )
    #     refseq_transcript_grch38 = relationship(
    #         "RefseqTranscripts",
    #         back_populates="refseq_transcripts_grch37",
    #     )
    #     # One to many relationships
    #     exons = relationship(
    #         "RefseqExons_Grch37", 
    #         back_populates="transcript"  # name of the exons column
    #     )
    #     cds = relationship(
    #         "RefseqCds_Grch37",
    #         back_populates="transcript"  # name of cds column
    #     )
    #     gene = relationship(
    #         "RefseqGenes_Grch37",
    #         back_populates="transcripts"
    #     )

    #     # Many to many relationships
    #     gencode_same_trans_coords_Grch37 = relationship(
    #         "Transcripts",
    #         secondary=association_table_trans_coords_Grch37,
    #         backref="gencode_transcripts"
    #     )
    #     gencode_same_exons_coords_Grch37 = relationship(
    #         "Transcripts",
    #         secondary=association_table_exon_coords_Grch37,
    #         backref="gencode_trans_same_exons"
    #     )
    #     gencode_same_cds_coords_Grch37 = relationship(
    #         "Transcripts",
    #         secondary=association_table_cds_coords_Grch37,
    #         backref="gencode_trans_same_cds"
    #     )
        
    #     def __repr__(self):
    #         return f"Refseq Transcript{self.id}"

    # class RefseqExons_Grch37(Base):
    #     __tablename__ = "RefseqExons_Grch37"
    #     exon_id = Column(String(50), primary_key=True)
    #     start = Column(Integer)
    #     end = Column(Integer)
    #     chromosome = Column(String(50))
    #     exon_number = Column(Integer)
    #     # Foreign key to relate with Transcript
    #     transcript_id = Column(
    #         String(50),
    #         ForeignKey("RefseqTranscripts_Grch37.id"),
    #         nullable=False
    #     )
    #     refseq_exon_grch38_id = Column(String(50), ForeignKey("RefseqExons.exon_id"))
    #     # Relationships
    #     refseq_exon_grch38 = relationship(
    #         "RefseqExons",
    #         back_populates="refseq_exon_grch37",
    #     )
    #     transcript = relationship(
    #         "RefseqTranscripts_Grch37",
    #         back_populates="exons"
    #     )

    #     def __repr__(self):
    #         return f"Refseq Exon {self.exon_id}"

    # class RefseqCds_Grch37(Base):
    #     __tablename__ = "RefseqCds_Grch37"
    #     primary_key = Column(Integer, primary_key=True, autoincrement=True)
    #     id = Column(String(50))
    #     start = Column(Integer)
    #     end = Column(Integer)
    #     chromosome = Column(String(50))
    #     # Foreign key to relate with Transcript
    #     transcript_id = Column(
    #         String(50),
    #         ForeignKey("RefseqTranscripts_Grch37.id")
    #     )

    #     transcript = relationship(
    #         "RefseqTranscripts_Grch37",
    #         back_populates="cds"
    #     )

    #     def __repr__(self):
    #         return f"Refseq Cds {self.id}"

    # class Transcripts_Grch37(Base):
    #     __tablename__ = "Transcripts_Grch37"
    #     transcript_id = Column(String(50), primary_key=True)
    #     mane_select = Column(Boolean, nullable=True, default=None)
    #     mane_clin = Column(Boolean, nullable=True, default=None)
    #     start = Column(Integer, nullable=False)
    #     end = Column(Integer, nullable=False)
    #     chromosome = Column(String(10), nullable=False)
    #     transcript_type = Column(String(50))
    #     transcript_name = Column(String(50))
    #     id = Column(String(50), nullable=False)
    #     havana_transcript = Column(String(50))
    #     protein_id = Column(String(50))
    #     gene_id = Column(String(50), ForeignKey("Genes.id"))
    #     lrg_transcript = Column(String(50))
    #     lrg_id = Column(String(50))
    #     ccds = Column(String(50))
    #     numb_exons = Column(Integer)
    #     numb_cds = Column(Integer)
    #     transcript_grch38_id = Column(
    #         String(50),
    #         ForeignKey("Transcripts.transcript_id")
    #     )

    #     transcript_grch38 = relationship(
    #         "Transcripts",
    #         back_populates="transcript_grch37",
    #     )

    #     # One to many relationships
    #     gene = relationship(
    #         "Genes_Grch37",
    #         back_populates="transcripts"  # column name of Genes relationship
    #     )
    #     exons = relationship(
    #         "Exons_Grch37",
    #         back_populates="transcript"  # column name of exons relationship
    #         )
    #     cds = relationship(
    #         "Cds_Grch37",
    #         back_populates="transcript"  # column name  of cds relationship
    #     )
    #     # Many to many relationships
    #     refseq_same_trans_coords_Grch37 = relationship(
    #         "RefseqTranscripts",
    #         secondary=association_table_trans_coords_Grch37,
    #         backref="refseq_transcripts_same_trans_coords",
    #     )
    #     refseq_same_exons_coords_Grch37 = relationship(
    #         "RefseqTranscripts",
    #         secondary=association_table_exon_coords_Grch37,
    #         backref="refseq_trans_same_exons_coords"
    #     )
    #     refseq_same_cds_coords_Grch37 = relationship(
    #         "RefseqTranscripts",
    #         secondary=association_table_cds_coords_Grch37,
    #         backref="refseq_trans_same_cds_coords"
    #     )

    #     def __repr__(self):
    #         return f"Transcript {self.transcript_id}"

    # class Exons_Grch37(Base):
    #     __tablename__ = "Exons_Grch37"
    #     exon_id = Column(String(50), primary_key=True)
    #     id = Column(String(50), nullable=False)
    #     exon_number = Column(Integer)
    #     start = Column(Integer, nullable=False)
    #     end = Column(Integer, nullable=False)
    #     chromosome = Column(String(50), nullable=False)
    #     transcript_id = Column(
    #         String(50),
    #         ForeignKey("Transcripts_Grch37.transcript_id"),
    #         nullable=True
    #     )
    #     exon_grch38_id = Column(String(50), ForeignKey("Exons.exon_id"))

    #     # Relationships
    #     exon_grch38 = relationship(
    #         "Exons",
    #         back_populates="exon_grch37",
    #     )

    #     transcript = relationship("Transcripts_Grch37", back_populates="exons")

    #     def __repr__(self):
    #         return f"Exon {self.exon_id}"

    # class Cds_Grch37(Base):
    #     __tablename__ = "Cds_Grch37"
    #     primary_key = Column(Integer, primary_key=True, autoincrement=True)
    #     id = Column(String(50))
    #     start = Column(Integer, nullable=False)
    #     end = Column(Integer, nullable=False)
    #     chromosome = Column(String(50), nullable=False)
    #     transcript_id = Column(
    #         String(50),
    #         ForeignKey("Transcripts_Grch37.transcript_id")
    #     )
    #     exon_id = Column(String(50))
    #     exon_number = Column(Integer)
    #     # Relationships
    #     transcript = relationship(
    #         "Transcripts_Grch37",
    #         back_populates="cds"
    #     )

    #     def __repr__(self):
    #         return f"CDS: {self.id}"

    # # inspector = inspect(engine)
    # # # Check if the Genes table exists in the database
    # # if not inspector.has_table("Genes"):
    # #     # Create all tables defined in the base class
    Base.metadata.create_all(engine)

    return (
        Gene_synonyms,
        Refseq_genes,
        Refseq_transcripts,
        Refseq_exons,
        Refseq_cds,
        Gencode_genes,
        Gencode_transcripts,
        Gencode_exons,
        Gencode_cds,
        File_version,
        refseq_gene_synonyms_association,
        gencode_gene_synonyms_association,
        association_table_trans_coords,
        association_table_exon_coords,
        association_table_cds_coords,
        session
    )


if "__main__" == __name__:
    (
        Gene_Synonyms,
        RefSeq_Genes,
        RefseqTranscripts,
        RefseqExons,
        RefseqCds,
        Genes,
        Transcripts,
        Exons,
        Cds,
        refseq_gene_synonyms_association,
        gencode_gene_synonyms_association,
        association_table_trans_coords,
        association_table_exons,
        association_table_cds,
        session
    ) = create_mysql_database()
