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
from typing import Tuple

logger = logging.getLogger(__name__)

class MissingGenomeFeatureFile3(Exception):
	pass

class InvalidGenome(Exception):
	pass

class GencodeRelease():
    '''
        Class that handles gencode human GFF data

        :param str genome_version: human genome version
        :param str output_dir: output directory
        :py:attr:`current_version`
        :py:meth:`download_gff3`
    '''
    def __init__(self, genome_version:str, output_dir:str):

        # Input validation
        valid_genomes = ['GRCh37', 'GRCh38']
        if genome_version not in valid_genomes:
            msg = (" ERROR: Invalid genome version {}. Choose between {}")\
                .format(genome_version, ','.join(valid_genomes))
            logging.error(msg)
            raise InvalidGenome(msg)

        self._genome_version = genome_version
        self._output_dir     = output_dir

        # Connect to ftp on instantiation
        self._origin_url = "ftp.ebi.ac.uk"
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

        self._working_dir    = os.path.join(self._output_dir, "GENCODE",
            self.current_release, self._genome_version)

        # Create working directory where all files will be stored
        Path(self._working_dir).mkdir(parents=True, exist_ok=True)

    @property
    def current_release(self) -> str:
        '''
            Returns the latest available release
        '''
		# move to pub dir and retrieve current readme
        self._ftp_cnn.cwd(self._parent_url)

        url = "pub/databases/gencode/Gencode_human/latest_release"
        self._ftp_cnn.cwd(url)
        files = self._ftp_cnn.nlst()

        target = next(filter(lambda x: x.startswith('gencode.v'), files), None)
        if not target:
            msg = " ERROR: missing file with gencode version number (gencode.v[0-9]+)"
            logging.error(msg)
            raise FileNotFoundError(msg)

        tmp_target = target.split(".")
        current_version = tmp_target[1].replace("v", "")
        print(current_version)

        return current_version

    def download_gff3(self) -> str:
        '''
            Download gff3 files from a genome version
            :param str output_dir: output directory. Not used
        '''
        # move back to parent url ("/")
        self._ftp_cnn.cwd(self._parent_url)

        current_v = self.current_release

        url_dict = {
            'GRCh37' : 'pub/databases/gencode/Gencode_human/latest_release/GRCh37_mapping/',
            'GRCh38' : 'pub/databases/gencode/Gencode_human/latest_release/'
        }

        gff_dict =  {
            'GRCh37' : f'gencode.v{current_v}lift37.annotation.gff3.gz',
            'GRCh38' : f'gencode.v{current_v}.annotation.gff3.gz'
        }

        print(self._genome_version)
        # move to genome version location
        self._ftp_cnn.cwd(url_dict[self._genome_version])
        dest_file = os.path.join(self._working_dir, gff_dict[self._genome_version])
        if not os.path.isfile(dest_file):
            msg = ("INFO: Downloading {}").format(gff_dict[self._genome_version])
            logging.info(msg)
            self._ftp_cnn.retrbinary("RETR " + gff_dict[self._genome_version], open(dest_file, 'wb').write)
        else:
            msg = ("INFO: File {} is already available at {}").format(gff_dict[self._genome_version], self._working_dir)
            logging.info(msg)
        self._ftp_cnn.quit()

        return dest_file
