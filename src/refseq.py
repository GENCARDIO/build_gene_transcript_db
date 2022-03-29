import os
import sys
import re
from pathlib import Path
from ftplib import FTP
import gzip
import tempfile
import logging
logger = logging.getLogger(__name__)


class MissingGenomeFeatureFile3(Exception):
	pass

class InvalidGenome(Exception):
	pass

refseq_to_chr = {
    "NC_000001" : "1",
    "NC_000002" : "2",
    "NC_000003" : "3",
    "NC_000004" : "4",
    "NC_000005" : "5",
    "NC_000006" : "6",
    "NC_000007" : "7",
    "NC_000008" : "8",
    "NC_000009" : "9",
    "NC_000010" : "0",
    "NC_000011" : "11",
    "NC_000012" : "12",
    "NC_000013" : "13",
    "NC_000014" : "14",
    "NC_000015" : "15",
    "NC_000016" : "16",
    "NC_000017" : "17",
    "NC_000018" : "18",
    "NC_000019" : "19",
    "NC_000020" : "20",
    "NC_000021" : "21",
    "NC_000022" : "22",
    "NC_000023" : "X",
    "NC_000024" : "Y"
}

class RefSeq():
	'''
		Class that handles RefSeq human GFF data
	'''
	def __init__(self, genome_version):
		self._genome_version = genome_version

		# Check that input genome version is valid
		valid_genomes = ['GRCh37', 'GRCh38', 'hg19', 'hg38']
		if not self._genome_version in valid_genomes:
			msg = "{} is not a valid genome. Accepted genomes are: {}"\
				.format(self._genome_version, ', '.join(valid_genomes))
			logging.error(msg)
			raise InvalidGenome(msg)

		# rename ucsc conventions to ensembl (hg* prefix to GRCh*)
		if self._genome_version == "hg19":
			self._genome_version = "GRCh37"
		if self._genome_version == "hg38":
			self._genome_version = "GRCh38"

		# Refseq base FTP
		self._refseq = "ftp.ncbi.nlm.nih.gov"

	@property
	def current_release(self):
		'''
		'''
		ftp = FTP("ftp.ncbi.nlm.nih.gov")
		ftp.login()
		latest_folder = "refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/"
		ftp.cwd(latest_folder)
		files = ftp.nlst()

		ori_file  = "README.txt"
		dest_file = ori_file

		if not os.path.isfile(dest_file):
			ftp.retrbinary("RETR " + ori_file, open(dest_file, 'wb').write)

		latest_refseq_human_release = ""
		with open (dest_file) as f:
			for line in f:
				line = line.rstrip("\n")
				tmp = line.split(" ")
				if line.startswith("name: Homo sapiens"):
					latest_refseq_human_release = tmp[-1]
		f.close()
		os.remove(dest_file)
		return latest_refseq_human_release

	def download_gff(self, output_dir):
		'''
		'''

		ftp = FTP(self._refseq)
		ftp.login()

		latest_folder = "refseq/H_sapiens/annotation/{}_latest/refseq_identifiers/"\
			.format(self._genome_version)
		ftp.cwd(latest_folder)

		ori_file  = "{}_latest_genomic.gff.gz".format(self._genome_version)
		dest_file = str(Path(output_dir)/ori_file)

		if not os.path.isfile(dest_file):
			msg = ("INFO: Downloading {}").format(ori_file)
			logging.info(msg)
			ftp.retrbinary("RETR " + ori_file, open(dest_file, 'wb').write)
		else:
			msg = ("INFO: File {} is already available at {}").format(ori_file, output_dir)
			logging.info(msg)
			ftp.quit()
			return dest_file
		ftp.quit()

		tmp_gff = dest_file.replace(".gff.gz", ".tmp.gff")
		o = open(tmp_gff, 'w')
		with gzip.open(dest_file, 'rt') as fin:
			for line in fin:
				line = line.rstrip("\n")
				if line.startswith("#"):
					o.write(line+"\n")
					continue
				tmp = line.split("\t")
				chr = tmp[0]
				tmp_chr = chr.split(".")
				base_chr = tmp_chr[0]
				if not base_chr in refseq_to_chr:
					continue
				tmp[0] = refseq_to_chr[base_chr]
				line = '\t'.join(tmp)
				o.write(line+"\n")
		o.close()

		tmp_gff_gz = ("{}.gz").format(tmp_gff)
		with open(tmp_gff, 'rb') as orig_file:
			with gzip.open(tmp_gff_gz, 'wb') as zipped_file:
				zipped_file.writelines(orig_file)
		os.remove(dest_file)
		os.remove(tmp_gff)
		os.rename(tmp_gff_gz, dest_file)
		return dest_file
