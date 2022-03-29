import os
import sys
import argparse
import logging
import re
from pathlib import Path
import logging
from datetime import datetime
from src.refseq import RefSeq
from src.ensembl import EnsemblRelease
from src.mane import Mane
from src.lrg import LRG
from src.database import build_database

def parse_args():

	parser = argparse.ArgumentParser(description="Build a gene/isoform database using ensembl and refseq gene models")
	parser.add_argument("--output_dir", help="Output directory", type=str, required=True)

	return parser.parse_args()

if __name__== "__main__":

	args = parse_args()
	output_dir     = args.output_dir
	Path(output_dir).mkdir(parents=True, exist_ok=True)

	now = datetime.now() # current date and time
	date_time = now.strftime("%d%m%Y")

	logging_file = output_dir + "/" + date_time + ".log"
	logging.basicConfig(format='%(asctime)s - %(message)s',
		filename=logging_file)
	logging.getLogger().setLevel(logging.INFO)
	logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

	genome_versions = ['GRCh38']
	versions_dict = {
		'ensembl': '.',
		'refseq' : '.',
		'mane'   : '.',
		'lrg'    : '.',
	}

	for genome_version in genome_versions:

		# MANE
		logging.info("INFO: Setting MANE transcripts")
		mane_main_dir = "MANE"
		Path(mane_main_dir).mkdir(parents=True, exist_ok=True)
		mane = Mane(genome_version)
		versions_dict['mane'] = mane.current_release
		mane_release_dir = mane_main_dir +"/"+versions_dict['mane'] + "/" + genome_version
		Path(mane_release_dir).mkdir(parents=True, exist_ok=True)
		mane_gff3 = mane.download_gff(mane_main_dir)

		# Ensembl
		logging.info("INFO: Setting Ensembl transcripts")
		ensembl_main_dir = "ENSEMBL"
		Path(ensembl_main_dir).mkdir(parents=True, exist_ok=True)
		ensembl = EnsemblRelease(genome_version)
		versions_dict['ensembl'] = ensembl.current_release
		ensembl_release_dir = ensembl_main_dir + "/" + versions_dict['ensembl']  +  "/" + genome_version
		Path(ensembl_release_dir).mkdir(parents=True, exist_ok=True)
		ensembl_gff3 = ensembl.download_gff(ensembl_release_dir)

		# RefSeq
		logging.info("INFO: Setting RefSeq transcripts")
		refseq_main_dir = "RefSeq"
		Path(refseq_main_dir).mkdir(parents=True, exist_ok=True)
		refseq = RefSeq(genome_version)
		versions_dict['refseq'] = refseq.current_release
		refseq_release_dir = refseq_main_dir + "/" + versions_dict['refseq'] +  "/" + genome_version
		Path(refseq_release_dir).mkdir(parents=True, exist_ok=True)
		refseq_gff3 = refseq.download_gff(refseq_release_dir)

	# LRG
	logging.info("INFO: Setting LRG transcripts")
	lrg_main_dir = "LRG"
	Path(lrg_main_dir).mkdir(parents=True, exist_ok=True)
	lrg = LRG()
	versions_dict['lrg'] = lrg.current_release
	lrg_release_dir = lrg_main_dir + "/" + versions_dict['lrg']
	Path(lrg_release_dir).mkdir(parents=True, exist_ok=True)
	lrg_file = lrg.download_transcripts(lrg_release_dir)

	input_file_dict = {
		'ensembl': ensembl_gff3,
		'refseq' : refseq_gff3,
		'mane'   : mane_gff3,
		'lrg'    : lrg_file
	}
	logging.info("INFO: Building gene-transcript database")
	build_database(input_file_dict, output_dir, versions_dict)
