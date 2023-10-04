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

from src.db import Base

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
