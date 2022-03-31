import os
import sys
import re
from pathlib import Path
from ftplib import FTP
import subprocess
import gzip
import logging
import glob
from collections import defaultdict
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

	@property
	def current_release(self) -> str:
		'''
			Returns the latest available release
		'''
		# move to pub dir and retrieve current readme
		current_readme = "current_README"
		ftp = FTP("ftp.ensembl.org")
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
			First, download GRCh38 gff3 as it is actively mantained.
			GRCh37 is not maintained since version 87 (2017).
			We found some discrepancies when converting 38 to 37 that are not optimal.
			For this reason, we sill keep old 87 release for GRCh37.

		'''

		ftp = FTP("ftp.ensembl.org")
		ftp.login()

		if self._genome_version == "GRCh37":

			# Move to homo sapiens repository
			ftp.cwd("pub/current_gff3/homo_sapiens/")
			files = ftp.nlst()
			# Check for the GRCh38 gff
			r = re.compile("Homo_sapiens\.GRCh38.\d+\.gff3.gz")
			gff_list = list(filter(r.match, files))
			if not gff_list:
				msg = ("ERROR: missing gff file from current Ensembl release")
				logging.error(msg)
				raise FileNotFoundError(msg)
			ori_file	= gff_list[0]
			grch38_gff = str(Path(output_dir)/ori_file)

			if not os.path.isfile(grch38_gff):
				msg = ("INFO: Downloading {}").format(ori_file)
				logging.info(msg)
				ftp.retrbinary("RETR " + ori_file, open(grch38_gff, 'wb').write)
			else:
				msg = ("INFO: File {} is already available at {}").format(ori_file, output_dir)
				logging.info(msg)
			ftp.quit()

			old_grch37_gff = self.download_old_grch37(output_dir)

			# Now convert GRCh38 to GRCh37
			lifted_grch37_gff = grch38_gff.replace("GRCh38", "GRCh37").replace(".gff3.gz", ".lift38.gff3.gz")
			grch37_gff = grch38_gff.replace("GRCh38", "GRCh37")
			if not os.path.isfile(lifted_grch37_gff):
				self.convert_38_to_37(grch38_gff, lifted_grch37_gff)

			# Note that some genes may be missing (e.g: CYP2D6)
			# As an ad-hoc solution we will mind the gaps using the old 87 release
			# self.build_grch37(old_grch37_gff, lifted_grch37_gff)

			os.remove(grch38_gff)
			return lifted_grch37_gff

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

	def sort_gff3(self, input_gff3):
		'''
			Sort a gff3
		'''
		output_gff3 = input_gff3.replace(".gff3", ".sorted.gff3")

		cmd = "sort -k1,1V -k4,4n -k5,5rn -k3,3r {} > {}".format(input_gff3, output_gff3)

		p1 = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
				stderr=subprocess.PIPE)
		output = p1.stdout.decode('UTF-8')
		error	= p1.stderr.decode('UTF-8')
		return output_gff3

	def coordinate_dict(self, input_gff3):
		'''
			From a gzipped gff3,
			Create a dict with keys pointing to unique coordinates
		'''
		coordinates_dict = defaultdict(dict)
		with gzip.open(input_gff3,'rt') as f:
			for line in f:
				line = line.rstrip("\n")
				if line.startswith("#"):
					continue
				tmp = line.split("\t")
				info = tmp[8].split(";")
				feature = tmp[2]
				chr  = tmp[0]
				pos  = tmp[3]
				end  = tmp[4]
				coordinates = "{}:{}-{}".format(chr,pos,end)
				if not coordinates in coordinates_dict:
					coordinates_dict[coordinates] = True
		f.close()
		return coordinates_dict

	def build_grch37(self, old_grch37_gff, lifted_grch37_gff):
		'''
		'''
		old_coordinates = self.coordinate_dict(old_grch37_gff)

		lifted_coordinates = self.coordinate_dict(lifted_grch37_gff)

		recatch_coords = defaultdict(dict)
		for old_coord in old_coordinates:
			if not old_coord in lifted_coordinates:
				recatch_coords[old_coord] = True

		temporary_grch37_gff = lifted_grch37_gff.replace(".gff3.gz", ".tmp.gff3")
		new_grch37_gff = lifted_grch37_gff.replace(".lift38.gff3.gz", ".gff3")

		output_dir = os.path.dirname(lifted_grch37_gff)
		only_87_gff = str(Path(output_dir)/"only.release.87.features.gff3")
		o = open(temporary_grch37_gff, 'w')
		p = open(only_87_gff, 'w')

		with gzip.open(old_grch37_gff,'rt') as f:
			for line in f:
				line = line.rstrip("\n")
				if line.startswith("#"):
					continue
				tmp = line.split("\t")
				chr  = tmp[0]
				pos  = tmp[3]
				end  = tmp[4]
				coordinates = "{}:{}-{}".format(chr,pos,end)
				if coordinates in recatch_coords:
					o.write(line+"\n")
					p.write(line+"\n")
		f.close()
		p.close()
		with gzip.open(lifted_grch37_gff,'rt') as f:
			for line in f:
				line = line.rstrip("\n")
				if line.startswith("#"):
					o.write(line+"\n")
				tmp = line.split("\t")
				chr  = tmp[0]
				pos  = tmp[3]
				end  = tmp[4]
				coordinates = "{}:{}-{}".format(chr,pos,end)
				if not coordinates in recatch_coords:
					o.write(line+"\n")
		f.close()
		o.close()

		sorted_gff3 = self.sort_gff3(temporary_grch37_gff)
		sorted_gff3_gz = ("{}.gz").format(sorted_gff3)

		with open(sorted_gff3, 'rb') as orig_file:
			with gzip.open(sorted_gff3_gz, 'wb') as zipped_file:
				zipped_file.writelines(orig_file)

		os.remove(old_grch37_gff)
		os.remove(temporary_grch37_gff)
		os.rename(sorted_gff3_gz, new_grch37_gff)

	def download_old_grch37(self, output_dir):
		'''
			Download legacy GRCh37 87 release)
		'''

		ftp = FTP("ftp.ensembl.org")
		ftp.login()
		ftp.cwd("pub/grch37/current/gff3/homo_sapiens/")
		files = ftp.nlst()

		r = re.compile("Homo_sapiens\.GRCh37.\d+\.gff3.gz")
		gff_list = list(filter(r.match, files))
		if not gff_list:
			msg = ("ERROR: missing gff file from current Ensembl release")
			logging.error(msg)
			raise FileNotFoundError(msg)

		ori_file  = gff_list[0]
		dest_file = str(Path(output_dir)/ori_file)
		if not os.path.isfile(dest_file):
			msg = ("INFO: Downloading legacy GRCh37 release {}").format(ori_file)
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
		msg = "INFO: {}".format(cmd)
		logging.info(msg)

		p1 = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
				stderr=subprocess.PIPE)
		output = p1.stdout.decode('UTF-8')
		error	= p1.stderr.decode('UTF-8')

		with open(output_file, 'rb') as orig_file:
			with gzip.open(output_gff, 'wb') as zipped_file:
				zipped_file.writelines(orig_file)
		os.remove(output_file)
