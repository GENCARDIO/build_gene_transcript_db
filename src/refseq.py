import os
import sys
import re
from pathlib import Path
from ftplib import FTP
import gzip
import tempfile
import logging
from typing import Union
logger = logging.getLogger(__name__)


class MissingGenomeFeatureFile3(Exception):
	pass

class InvalidGenome(Exception):
	pass

class RefSeq():
	'''
		Class that handles RefSeq human GFF data
        :param str genome_version: human genome version
        :param str output_dir: output directory
        :py:attr:`current_version`
        :py:meth:`download_gff3`
	'''
	def __init__(self, genome_version:str, output_dir:str):

		valid_genomes = ['GRCh37', 'GRCh38']
		if genome_version not in valid_genomes:
			msg = (" ERROR: Invalid genome version {}. Choose between {}")\
				.format(genome_version, ','.join(valid_genomes))
			logging.error(msg)
			raise InvalidGenome(msg)
		self._genome_version = genome_version

		# Connect to ftp on instantiation
		self._origin_url = "ftp.ncbi.nlm.nih.gov"
		try:
			self._ftp_cnn = FTP(self._origin_url)
			self._ftp_cnn.login()
		except:
			msg = ("ERROR: Could not stablish a connection with {}")\
				.format(self._origin_url)
			logging.error(msg)
			raise ConnectionError(msg)
		else:
			msg = ("INFO: Connected to {}").format(self._origin_url )
			logging.info(msg)
		self._parent_url = self._ftp_cnn.pwd()

		self._output_dir     = output_dir
		self._working_dir    = os.path.join(self._output_dir, "RefSeq",
			self.current_release, self._genome_version)

		# Create working directory where all files will be stored
		Path(self._working_dir).mkdir(parents=True, exist_ok=True)

	@property
	def current_release(self) -> str:
		'''
		'''
		# move back to parent url ("/")
		self._ftp_cnn.cwd(self._parent_url)

		url = "refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers"
		try:
			self._ftp_cnn.cwd(url)
		except:
			msg = ("ERROR: Unable to move to {}").format(url)
			logging.error(msg)
			raise ConnectionError(msg)

		ori_file  = "README.txt"
		dest_file = ori_file

		if not os.path.isfile(dest_file):
			self._ftp_cnn.retrbinary("RETR " + ori_file, open(dest_file, 'wb').write)

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

	def download_gff3(self) -> str:
		'''
			Download gff3 files from a genome version
			:param str output_dir: output directory. Not used
		'''
		# move back to parent url ("/")
		self._ftp_cnn.cwd(self._parent_url)

		url = "refseq/H_sapiens/annotation/{}_latest/refseq_identifiers"\
			.format(self._genome_version)
		try:
			self._ftp_cnn.cwd(url)
		except:
			msg = ("ERROR: Unable to move to {}").format(url)
			logging.error(msg)
			raise ConnectionError(msg)

		ori_file  = ("{}_latest_genomic.gff.gz").format(self._genome_version)
		dest_file = os.path.join(self._working_dir, ori_file)

		if not os.path.isfile(dest_file):
			msg = ("INFO: Downloading {}").format(ori_file)
			logging.info(msg)
			self._ftp_cnn.retrbinary("RETR " + ori_file, open(dest_file, 'wb').write)
		else:
			msg = ("INFO: File {} is already available at {}").format(ori_file,
				self._working_dir,)
			logging.info(msg)
			self._ftp_cnn.quit()
			return dest_file
		self._ftp_cnn.quit()

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
				refseq_chr = tmp_chr[0]

				chrom = self.refseq_to_chr_convention(refseq_chr)
				if chrom is None:
					continue
				tmp[0] = chrom
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

	@staticmethod
	def refseq_to_chr_convention(self, chrom: str) -> Union[str, None]:
		'''
			Convert a refseq chromosome sequence to its molecule name homologue
			:param str chrom: input refesq chromosome
		'''
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
		if not chrom in refseq_to_chr:
			return None

		return refseq_to_chr[chrom]
