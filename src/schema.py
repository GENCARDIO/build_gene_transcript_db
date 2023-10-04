from pydantic import BaseModel

class Gencode_genes(BaseModel):
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