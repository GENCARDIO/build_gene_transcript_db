from dataclasses import dataclass

@dataclass
class EnsemblGenes:
    gene_id: str
    start: int
    end: int
    chromosome: str
    id: str
    tag: str
    num_gencode_transcripts: int
    num_refseq_transcripts: int
    refseq_same_gene_coords: str
    gene_name: str
    gene_type: str
    level: int
    havana_gene: str
    hgnc_id: int

@dataclass
class RefseqGene:
    gene_id: str
    start: int
    end: int
    chromosome: int
    id: str
    hgnc_id: int
    name: str
    gene_biotype: str

@dataclass
class EnsemblTranscript:
    


