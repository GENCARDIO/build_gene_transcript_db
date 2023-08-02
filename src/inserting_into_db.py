from src.global_variables import logging
from sqlalchemy.exc import SQLAlchemyError

batch_size = 1000

def insert_db_gencode_gene_data(gencode_genename_geneobject, session, Genes, Gene_Synonyms, gencode_gene_synonyms_association, RefSeq_Genes):
    """
    From a dictionary with [gene_id]:[feature]:[values]
    inserts the features and values into the genes table

    Params:
        gencode_gene_dict[gene_id]:[feature]:[values]
        session: session of the database, returned by function create_database
    """
    genes_inserted = 0
    for gene_name, gene_items in gencode_genename_geneobject.items():
        for ensembl_gene in gene_items:
            synonyms = []
            refseq_same_gene_coords = []
            for gene_synonym in ensembl_gene.gene_synonyms:
                existing_synonym = session.query(Gene_Synonyms).filter_by(synonym=gene_synonym).first()
                if existing_synonym is None:
                    new_synonym = Gene_Synonyms(synonym=gene_synonym)
                    session.add(new_synonym)
                    synonyms.append(new_synonym)
                else:
                    synonyms.append(existing_synonym)
                session.commit()
            if ensembl_gene.refseq_same_gene_coords is not None:
                for refseq_same_gene_coord in ensembl_gene.refseq_same_gene_coords:
                    refseq_db_ins = session.query(RefSeq_Genes).filter_by(primary_key=refseq_same_gene_coord.primary_key).first()
                    refseq_same_gene_coords.append(refseq_db_ins)
                existing_synonym = session.query(Gene_Synonyms).filter_by(synonym=gene_synonym).first()
                ensembl_gene.refseq_same_gene_coords = list(ensembl_gene.refseq_same_gene_coords)
            else:
                ensembl_gene.refseq_same_gene_coords = []
            logging.info(f"inserted {ensembl_gene}")
            print(ensembl_gene.refseq_same_gene_coords)
            new_gene = Genes(
                gene_id_genome_version=ensembl_gene.get_full_id(),
                gene_id=ensembl_gene.gene_id,
                genome_version=ensembl_gene.genome_version,
                start=ensembl_gene.start,
                end=ensembl_gene.end,
                chromosome=ensembl_gene.chromosome,
                id=ensembl_gene.id,
                gene_name=ensembl_gene.gene_name,
                gene_type=ensembl_gene.gene_type,
                level=ensembl_gene.level,
                synonyms=synonyms,
                tag=ensembl_gene.tag,
                havana_gene=ensembl_gene.havana_gene,
                hgnc_id=ensembl_gene.hgnc_id,
                num_gencode_transcripts=ensembl_gene.num_gencode_transcripts,
                liftover_chromosome=ensembl_gene.liftover_chromosome,
                liftover_start=ensembl_gene.liftover_start,
                liftover_end=ensembl_gene.liftover_end
            )

        # new_refseq_gene = RefSeq_Genes(
        #     gene_id=refseq_gene.gene_id,
        #     start=refseq_gene.start,
        #     end=refseq_gene.end,
        #     chromosome=refseq_gene.chromosome,
        #     id=refseq_gene.id,
        #     name=refseq_gene.name,
        #     gene_biotype=refseq_gene.gene_biotype,
        #     feature=refseq_gene.feature,
        #     gene_synonyms=refseq_gene.gene_synonyms,
        #     hgnc_id=refseq_gene.hgnc_id,
        #     num_transcripts=refseq_gene.num_transcripts,
        #     gencode_same_gene_coords=refseq_gene.gencode_same_gene_coords

        # )
        try:
            print(f"inserting gene {new_gene}")
            # Add new gene to the session
            session.add(new_gene)
            genes_inserted += 1
            if genes_inserted % batch_size == 0:
                session.commit()
                logging.info(f"{genes_inserted} genes inserted into db")
        except:
            session.rollback()
            msg = "Error in database writing"
            raise SQLAlchemyError(msg)
        else:
            msg = f"gene {new_gene} record successfully"
            logging.info(msg)

    if genes_inserted % batch_size != 0:
        session.commit()
        logging.info(f"{genes_inserted} genes into db")


def insert_db_refseq_genes(refseq_genename_geneobject, session, RefSeq_Genes, Gene_Synonyms):
    gene_inserted = 0
    for gene_name, gene_items in refseq_genename_geneobject.items():
        for refseq_gene in gene_items:
            synonyms = []
            new_synonym = None
            if refseq_gene.gene_synonyms is not None:
                for gene_synonym in refseq_gene.gene_synonyms:
                    existing_synonym = session.query(Gene_Synonyms).filter_by(synonym=gene_synonym).first()
                    if existing_synonym is None:
                        new_synonym = Gene_Synonyms(synonym=gene_synonym)
                        synonyms.append(new_synonym)
                        session.add(new_synonym)
                    else:
                        synonyms.append(existing_synonym)
                    
                session.commit()
            if refseq_gene.gencode_same_gene_coords is not None:
                refseq_gene.gencode_same_gene_coords = list(refseq_gene.gencode_same_gene_coords)
            else:
                refseq_gene.gencode_same_gene_coords = []
            logging.info(f"inserted {refseq_gene}")
            print(refseq_gene.gencode_same_gene_coords)
            new_gene = RefSeq_Genes(
                primary_key=refseq_gene.primary_key,
                gene_id=refseq_gene.gene_id,
                genome_version=refseq_gene.genome_version,
                start=refseq_gene.start,
                end=refseq_gene.end,
                chromosome=refseq_gene.chromosome,
                gene_name=refseq_gene.name,
                gene_type=refseq_gene.gene_type,
                synonyms=synonyms,
                hgnc_id=refseq_gene.hgnc_id,
                num_transcripts=refseq_gene.num_transcripts,
                liftover_chromosome=refseq_gene.liftover_chromosome,
                liftover_start=refseq_gene.liftover_start,
                liftover_end=refseq_gene.liftover_end
            )
            try:
                print(f"inserting gene {new_gene}")
                # Add new gene to the session
                session.add(new_gene)
                gene_inserted += 1
                if gene_inserted % batch_size == 0:
                    session.commit()
                    logging.info(f"{gene_inserted} genes into db")

            except:
                session.rollback()
                msg = "Error in database writing"
                raise SQLAlchemyError(msg)
            else:
                msg = f"gene {new_gene} record successfully"
                logging.info(msg)
            
            refseq_gene.primary_key = new_gene.primary_key

    if gene_inserted % batch_size != 0:
        session.commit()
        logging.info(f"{gene_inserted} genes into db")
def insert_db_refseq_transcripts(refseq_genename_geneid_transobject, session, RefseqTranscripts):
    transcripts_inserted = 0
    for gene_name, genes_dict in refseq_genename_geneid_transobject.items():
        for gene_obj, transobjects in genes_dict.items():
            for transobj in transobjects:
                new_transcript = RefseqTranscripts(
                    id_genome_version=transobj.get_full_id(),
                    id=transobj.id,
                    genome_version=transobj.genome_version,
                    start=transobj.start,
                    end=transobj.end,
                    chromosome=transobj.chromosome,
                    gene_id=transobj.gene_id,
                    gene_name=transobj.gene_name,
                    hgnc_id=transobj.hgnc_id,
                    lrg_id=transobj.lrg_id,
                    lrg_transcript=transobj.lrg_transcript,
                    ccds=transobj.ccds,
                    numb_exons=transobj.numb_exons,
                    numb_cds=transobj.numb_cds,
                    mane_clin=transobj.mane_clin,
                    mane_select=transobj.mane_select,
                    gene_pk=gene_obj.primary_key,
                    liftover_chromosome=transobj.liftover_chromosome,
                    liftover_start=transobj.liftover_start,
                    liftover_end=transobj.liftover_end
                )
                try:
                    print(f"inserting gene {new_transcript}")
                    # Add new gene to the session
                    session.add(new_transcript)
                    transcripts_inserted += 1
                    if transcripts_inserted % batch_size == 0:
                        session.commit()
                        logging.info(f"{transcripts_inserted} transcripts into db")
                except:
                    session.rollback()
                    msg = "Error in database writing"
                    raise SQLAlchemyError(msg)
                else:
                    msg = f"gene {new_transcript} record successfully"
                    logging.info(msg)
    if transcripts_inserted % batch_size != 0:
        session.commit()
        logging.info(f"{transcripts_inserted} transcripts into db")

def insert_db_refseq_exons(refseq_genename_geneid_transid_exonobject, session, RefseqExons):
    exons_inserted = 0
    for gene_name, genes_dict in refseq_genename_geneid_transid_exonobject.items():
        for gene_id, trans_dict in genes_dict.items():
            for trans_obj, exon_objects in trans_dict.items():
                for exon_obj in exon_objects:
                    new_exon = RefseqExons(
                        exon_id_genome_version=exon_obj.get_full_id(),
                        exon_id=exon_obj.exon_id,
                        genome_version=exon_obj.genome_version,
                        start=exon_obj.start,
                        end=exon_obj.end,
                        chromosome=exon_obj.chromosome,
                        exon_number=exon_obj.exon_number,
                        transcript_id=exon_obj.transcript_id,
                        liftover_chromosome=exon_obj.liftover_chromosome,
                        liftover_start=exon_obj.liftover_start,
                        liftover_end=exon_obj.liftover_end
                    )
                    print(exon_obj.transcript_id, "transcript_id")
                    try:
                        print(f"inserting gene {new_exon}")
                        # Add new gene to the session
                        session.add(new_exon)
                        exons_inserted += 1
                        if exons_inserted % batch_size == 0:
                            session.commit()
                            logging.info(f"{exons_inserted} exons into db")
                    except:
                        session.rollback()
                        msg = "Error in database writing"
                        raise SQLAlchemyError(msg)
                    else:
                        msg = f"gene {new_exon} record successfully"
                        logging.info(msg)
    if exons_inserted % batch_size != 0:
        session.commit()
        logging.info(f"{exons_inserted} exons into db")
def insert_db_refseq_cds(refseq_genename_geneid_transid_cdsobject, session, RefseqCds):
    cds_inserted = 0
    for gene_name, genes_dict in refseq_genename_geneid_transid_cdsobject.items():
        for gene_id, trans_dict in genes_dict.items():
            for trans_obj, cds_objects in trans_dict.items():
                for cds_obj in cds_objects:
                    if trans_obj.id != cds_obj.transcript_id:
                        raise(ValueError(f"transcript id from transcript object and transcript id from cds object don't match, transcript {trans_obj.id} cds: {cds_obj.transcript_id}"))

                    new_cds = RefseqCds(
                        primary_key=cds_obj.primary_key,
                        id=cds_obj.id,
                        genome_version=cds_obj.genome_version,
                        start=cds_obj.start,
                        end=cds_obj.end,
                        chromosome=cds_obj.chromosome,
                        transcript_id=cds_obj.transcript_id,
                        liftover_chromosome=cds_obj.liftover_chromosome,
                        liftover_start=cds_obj.liftover_start,
                        liftover_end=cds_obj.liftover_end
                    )
                    print(cds_obj.transcript_id, "transcript_id")
                    try:
                        print(f"inserting gene {new_cds}")
                        # Add new gene to the session
                        session.add(new_cds)
                        cds_inserted += 1
                        if cds_inserted % batch_size == 0:
                            session.commit()
                            logging.info(f"{cds_inserted} cds into db")
                    except:
                        session.rollback()
                        msg = "Error in database writing"
                        raise SQLAlchemyError(msg)
                    else:
                        msg = f"gene {new_cds} record successfully"
                        logging.info(msg)
        
        if cds_inserted % batch_size != 0:
            session.commit()
            logging.info(f"{cds_inserted} cds into db")


def insert_db_gencode_transcripts(gencode_genename_geneid_transobject, session, Transcripts, RefseqTranscripts):
    transcript_inserted = 0
    for gene_name, genes_dict in gencode_genename_geneid_transobject.items():
        for gene_obj, transobjects in genes_dict.items():
            for transobj in transobjects:
                refseq_db_same_trans_coords = []
                refseq_db_same_exons_coords = []
                refseq_db_same_cds_coords = []
                print(transobj.refseq_same_trans_coords, transobj.refseq_same_exons_coords, transobj.refseq_same_cds_coords)
                if transobj.refseq_same_trans_coords:
                    for refseq_id in transobj.refseq_same_trans_coords:
                        refseq_db_ins = session.query(RefseqTranscripts).filter_by(id_genome_version=refseq_id).first()
                        refseq_db_same_trans_coords.append(refseq_db_ins)
                if transobj.refseq_same_exons_coords:
                    for refseq_id in transobj.refseq_same_exons_coords:
                        refseq_db_trans_ins = session.query(RefseqTranscripts).filter_by(id_genome_version=refseq_id).first()
                        refseq_db_same_exons_coords.append(refseq_db_trans_ins)
                if transobj.refseq_same_cds_coords:
                    for refseq_id in transobj.refseq_same_cds_coords:
                        refseq_db_trans_cds_ins = session.query(RefseqTranscripts).filter_by(id_genome_version=refseq_id).first()
                        refseq_db_same_cds_coords.append(refseq_db_trans_cds_ins)
                new_transcript = Transcripts(
                    transcript_id_genome_version=transobj.get_full_id(),
                    transcript_id=transobj.transcript_id,
                    genome_version=transobj.genome_version,
                    mane_select=transobj.mane_select,
                    mane_clin=transobj.mane_clin,
                    start=transobj.start,
                    end=transobj.end,
                    chromosome=transobj.chromosome,
                    transcript_type=transobj.transcript_type,
                    transcript_name=transobj.transcript_name,
                    id=transobj.id,
                    havana_transcript=transobj.havana_transcript,
                    protein_id=transobj.protein_id,
                    gene_id=transobj.parent,
                    lrg_transcript=transobj.lrg_transcript,
                    lrg_id=transobj.lrg_id,
                    ccds=transobj.ccds,
                    numb_exons=transobj.numb_exons,
                    numb_cds=transobj.numb_cds,
                    liftover_chromosome=transobj.liftover_chromosome,
                    liftover_start=transobj.liftover_start,
                    liftover_end=transobj.liftover_end,
                    refseq_same_trans_coords=refseq_db_same_trans_coords,
                    refseq_same_exons_coords=refseq_db_same_exons_coords,
                    refseq_same_cds_coords=refseq_db_same_cds_coords
                )
                try:
                    print(f"inserting gene {new_transcript}")
                    # Add new gene to the session
                    session.add(new_transcript)
                    transcript_inserted += 1
                    if transcript_inserted % batch_size == 0:
                        session.commit()
                        logging.info(f"{transcript_inserted} transcripts into db")
                except:
                    session.rollback()
                    msg = "Error in database writing"
                    raise SQLAlchemyError(msg)
                else:
                    msg = f"gene {new_transcript} record successfully"
                    logging.info(msg)
    
    if transcript_inserted % batch_size != 0:
        session.commit()
        logging.info(f"{transcript_inserted} transcripts into db")
    
def insert_db_gencode_exons(gencode_genename_geneid_transid_exonobject, session, Exons):
    exons_inserted = 0
    for gene_name, genes_dict in gencode_genename_geneid_transid_exonobject.items():
        for gene_id, trans_dict in genes_dict.items():
            for trans_obj, exon_objects in trans_dict.items():
                for exon_obj in exon_objects:
                    new_exon = Exons(
                        exon_id_genome_version=exon_obj.get_full_id(),
                        exon_id=exon_obj.id,
                        genome_version=exon_obj.genome_version,
                        id=exon_obj.exon_id,
                        exon_number=exon_obj.exon_number,
                        start=exon_obj.start,
                        end=exon_obj.end,
                        chromosome=exon_obj.chromosome,
                        transcript_id=trans_obj.get_full_id(),
                        liftover_chromosome=exon_obj.liftover_chromosome,
                        liftover_start=exon_obj.liftover_start,
                        liftover_end=exon_obj.liftover_end
                    )
                    print(exon_obj.transcript_id, "transcript_id")
                    try:
                        print(f"inserting exon {new_exon}")
                        # Add new gene to the session
                        session.add(new_exon)
                        exons_inserted += 1
                        if exons_inserted % batch_size == 0:
                            session.commit()
                            logging.info(f"{exons_inserted} exons inserted correctly into db")
                    except:
                        session.rollback()
                        msg = "Error in database writing"
                        raise SQLAlchemyError(msg)
                    else:
                        msg = f"exon {new_exon} record successfully"
                        logging.info(msg)
        if exons_inserted % batch_size != 0:
            session.commit()
            logging.info(f"{exons_inserted} exons inserted successfully into db")


def insert_db_gencode_cds(gencode_genename_geneid_transid_cdsobject, session, Cds):
    cds_inserted = 0
    for gene_name, genes_dict in gencode_genename_geneid_transid_cdsobject.items():
        for gene_id, trans_dict in genes_dict.items():
            for trans_obj, cds_objects in trans_dict.items():
                for cds_obj in cds_objects:
                    # if trans_obj.id != cds_obj.parent:
                    #     raise(ValueError(f"transcript id from transcript object and transcript id from cds object don't match, transcript {trans_obj.id} cds: {cds_obj.transcript_id}"))
                    new_cds = Cds(
                        primary_key=cds_obj.primary_key,
                        id=cds_obj.id,
                        genome_version=cds_obj.genome_version,
                        start=cds_obj.start,
                        end=cds_obj.end,
                        chromosome=cds_obj.chromosome,
                        exon_id=cds_obj.exon_id,
                        exon_number=cds_obj.exon_number,
                        transcript_id=trans_obj.get_full_id(),
                        liftover_chromosome=cds_obj.liftover_chromosome,
                        liftover_start=cds_obj.liftover_start,
                        liftover_end=cds_obj.liftover_end
                    )
                    print(new_cds.transcript_id, trans_obj.get_full_id())
                    try:
                        print(f"inserting cds {new_cds}")
                        # Add new gene to the session
                        session.add(new_cds)
                        cds_inserted += 1
                        # commit changes in batches of 1000
                        if cds_inserted % batch_size == 0:
                            session.commit()
                            logging.info(f"{cds_inserted} cds inserted successfully")
                    except:
                        session.rollback()
                        msg = "Error in database writing"
                        raise SQLAlchemyError(msg)
                    else:
                        msg = f"cds {new_cds} record successfully"
                        logging.info(msg)
    # Commit any remaining changes that are not part of a complete batch
    if cds_inserted % batch_size != 0:
        session.commit()
        logging.info(f"{cds_inserted} genes inserted successfully")