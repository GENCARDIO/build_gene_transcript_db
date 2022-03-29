import os
import sys
import re
from pathlib import Path
from ftplib import FTP
import subprocess
import gzip
import logging
logger = logging.getLogger(__name__)


class MissingGenomeFeatureFile3(Exception):
		pass

class InvalidGenome(Exception):
		pass

class EnsemblRelease():
	'''
			Class that handles ensembl human GFF data
	'''

	def __init__(self, genome_version, release="current"):
		self._genome_version = genome_version
		self._current = release

		# Check that input genome version is valid
		valid_genomes = ['GRCh37', 'GRCh38', 'hg19', 'hg38']
		if not self._genome_version in valid_genomes:
			msg = "{} is not a valid genome. Accepted genomes are: {}"\
					.format(self._genome_version, ', '.join(valid_genomes))
			logging.error(msg)
			raise InvalidGenome(msg)

		# rename ucsc conventions to ensembl (hg prefix to GRCh)
		if self._genome_version == "hg19":
			self._genome_version = "GRCh37"
		if self._genome_version == "hg38":
			self._genome_version = "GRCh38"

		# Ensembl base FTP
		self._ensembl_ftp = "ftp.ensembl.org"

	@property
	def current_release(self) -> str:
		'''
			Returns the latest available release
		'''
		# move to pub dir and retrieve current readme
		current_readme = "current_README"
		ftp = FTP(self._ensembl_ftp)
		ftp.login()
		ftp.cwd("pub/")
		ftp.retrbinary("RETR " + current_readme, open(current_readme, 'wb').write)
		ftp.quit()

		release_version = ""
		with open(current_readme) as f:
			for line in f:
				line = line.rstrip("\n")

				# version when lines starts with Ensembl Release
				if line.startswith("Ensembl Release"):
					tmp = line.split(" ")
					release_version = tmp[2]
		f.close()
		if not release_version:
			msg = ("ERROR: Ensembl current release version was not found")
			logging.error(msg)
			sys.exit()
		os.remove(current_readme)

		return release_version

	def download_gff(self, output_dir) -> str:
		'''
			Automatically downloads GRCh38 gff3.
			For version GRCh37, it converts coordinates back from GRCh38
		'''

		ftp = FTP(self._ensembl_ftp)
		ftp.login()

		if self._genome_version == "GRCh37":

			ftp.cwd("pub/current_gff3/homo_sapiens/")
			files = ftp.nlst()

			r = re.compile("Homo_sapiens\.GRCh38.\d+\.gff3.gz")
			gff_list = list(filter(r.match, files))

			if not gff_list:
				msg = ("ERROR: missing gff file from current Ensembl release")
				logging.error(msg)
				raise FileNotFoundError(msg)

			ori_file	= gff_list[0]
			dest_file = str(Path(output_dir)/ori_file)
			output_gff = dest_file.replace("GRCh38", "GRCh37")
			if not os.path.isfile(dest_file):
				msg = ("INFO: Downloading {}").format(ori_file)
				logging.info(msg)
				ftp.retrbinary("RETR " + ori_file, open(dest_file, 'wb').write)
			else:
				msg = ("INFO: File {} is already available at {}").format(ori_file, output_dir)
				logging.info(msg)
			ftp.quit()

			output_gff = dest_file.replace("GRCh38", "GRCh37")
			if not os.path.isfile(output_gff):
				self.convert_38_to_37(dest_file, output_gff)
			os.remove(dest_file)

			return output_gff

		if self._genome_version == "GRCh38":

			ftp.cwd("pub/current_gff3/homo_sapiens/")
			files = ftp.nlst()

			r = re.compile("Homo_sapiens\.GRCh38.\d+\.gff3.gz")
			gff_list = list(filter(r.match, files))

			if not gff_list:
				msg = ("ERROR: missing gff file from current Ensembl release")
				logging.error(msg)
				raise FileNotFoundError(msg)

			ori_file	= gff_list[0]
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

	def convert_38_to_37(self, input_gff, output_gff):
		'''
		'''
		# Issue with hg19/GRCh37: this genome version is not actively being
		chain = "chains/hg38ToHg19.over.chain.gz"
		output_file = output_gff.replace(".gz", "")

		cmd = ('CrossMap.py {} {} {} {}').format("gff", chain, input_gff, output_file)
		msg = "INFO: " + cmd
		logging.info(msg)

		p1 = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
				stderr=subprocess.PIPE)
		output = p1.stdout.decode('UTF-8')
		error	= p1.stderr.decode('UTF-8')

		with open(output_file, 'rb') as orig_file:
				with gzip.open(output_gff, 'wb') as zipped_file:
						zipped_file.writelines(orig_file)
		os.remove(output_file)
