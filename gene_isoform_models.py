import os
import sys
import argparse
import logging
import re
from pathlib import Path
import logging
from collections import defaultdict
from datetime import datetime
from src.gencode import GencodeRelease
from src.refseq import RefSeq
# from src.ensembl import EnsemblRelease
from src.mane import Mane
from src.lrg import LRG
from src.database import build_database
from src.utilities import check_crossmap

def parse_args():
	'''
	'''
	parser = argparse.ArgumentParser(description="Build a gene/isoform database using ensembl/refseq gene models")
	parser.add_argument("--output_dir", help="Output directory", type=str, required=True)
	return parser.parse_args()

if __name__== "__main__":

	args = parse_args()
	output_dir     = args.output_dir
	Path(output_dir).mkdir(parents=True, exist_ok=True)

	now = datetime.now() # current date and time
	date_time = now.strftime("%d%m%Y") + ".log"

	logging_file = os.path.join(output_dir, date_time)
	logging.basicConfig(format='%(asctime)s - %(message)s',
		filename=logging_file)
	logging.getLogger().setLevel(logging.INFO)
	logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

	# genome_versions = ['GRCh37', 'GRCh38']
	genome_versions = ['GRCh37']
	versions_dict = defaultdict(dict)
	versions_dict = {
		'gencode': '.',
		'ensembl': '.',
		'refseq' : '.',
		'mane'   : '.',
		'lrg'    : '.',
	}
	gff_files_dict = {
		'gencode': '.',
		'ensembl': '.',
		'refseq' : '.',
		'mane'   : '.',
		'lrg'    : '.'
	}

	# Check that CrossMap is available throught PATH
	check_crossmap()

	# Now do the actual job
	for genome_version in genome_versions:

		versions_dict[genome_version]  = defaultdict(dict)
		gff_files_dict[genome_version] = defaultdict(dict)

		# GENCODE
		msg = ("INFO: Downloading GENCODE {} transcript annotation").format(genome_version)
		logging.info(msg)
		gencode = GencodeRelease(genome_version, output_dir)
		versions_dict[genome_version]['gencode']  = gencode.current_release
		gencode_gff3 = gencode.download_gff3()
		gff_files_dict[genome_version]['gencode'] = gencode_gff3

		# RefSeq
		msg = ("INFO: Downloading RefSeq transcripts for version {}").format(genome_version)
		logging.info(msg)
		refseq = RefSeq(genome_version, output_dir)
		versions_dict[genome_version]['refseq']  = refseq.current_release
		refseq_gff3 = refseq.download_gff3()
		gff_files_dict[genome_version]['refseq'] = refseq_gff3

		# LRG
		msg = ("INFO: Downloading LRG transcripts for version {}").format(genome_version)
		logging.info(msg)
		lrg = LRG(output_dir)
		versions_dict[genome_version]['lrg']  = lrg.current_release
		lrg_file = lrg.download_transcripts()
		gff_files_dict[genome_version]['lrg'] = lrg_file

		# MANE
		msg = ("INFO: Downloading MANE transcripts for version {}").format(genome_version)
		logging.info(msg)
		mane = Mane(genome_version, output_dir)
		versions_dict[genome_version]['mane']  = mane.current_release
		mane_gff3 = mane.download_gff3()
		gff_files_dict[genome_version]['mane'] = mane_gff3

	logging.info("INFO: Building gene-transcript database")
	build_database("GRCh37", gff_files_dict, output_dir, versions_dict)
