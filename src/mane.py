import os
import sys
import re
from pathlib import Path
from ftplib import FTP
import gzip
import subprocess
import logging
logger = logging.getLogger(__name__)


class Mane():
    '''
        Class that handles MANE human GFF data
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
        self._output_dir = output_dir
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

        self._working_dir    = os.path.join(self._output_dir, "MANE",
            self.current_release)
        Path(self._working_dir).mkdir(parents=True, exist_ok=True)
        self._curr_release = self.current_release

    @property
    def current_release(self)->str:
        '''
        '''
        self._ftp_cnn.cwd(self._parent_url)

        url = "refseq/MANE/MANE_human/current"
        try:
            self._ftp_cnn.cwd(url)
        except:
            msg = ("ERROR: Unable to move to {}").format(url)
            logging.error(msg)
            raise ConnectionError(msg)

        # latest_folder = "refseq/MANE/MANE_human/current/"
        # ftp.cwd(latest_folder)
        ori_file  = "README_versions.txt"
        dest_file = ori_file

        if not os.path.isfile(dest_file):
            self._ftp_cnn.retrbinary("RETR " + ori_file, open(dest_file, 'wb').write)

        latest_mane_version = ""
        with open (dest_file) as f:
            for line in f:
                line = line.rstrip("\n")
                tmp = line.split("\t")
                if line.startswith("MANE"):
                    latest_mane_version = tmp[1]
        f.close()
        os.remove(dest_file)
        return latest_mane_version

    def download_gff3(self) -> str:
        '''
        '''

        self._ftp_cnn.cwd(self._parent_url)

        url = "refseq/MANE/MANE_human/current"
        try:
            self._ftp_cnn.cwd(url)
        except:
            msg = ("ERROR: Unable to move to {}").format(url)
            logging.error(msg)
            raise ConnectionError(msg)

        grch38_dir =  os.path.join(self._working_dir, "GRCh38")
        Path(grch38_dir).mkdir(parents=True, exist_ok=True)

        grch38_gff_name  = "MANE.GRCh38.v{}.ensembl_genomic.gff.gz".format(self._curr_release)
        grch38_gff = str(Path(grch38_dir)/grch38_gff_name)
        if not os.path.isfile(grch38_gff):
            self._ftp_cnn.retrbinary("RETR " + grch38_gff_name, open(grch38_gff, 'wb').write)

        grch37_dir =  os.path.join(self._working_dir, "GRCh37")
        Path(grch37_dir).mkdir(parents=True, exist_ok=True)
        grch37_gff_name = grch38_gff_name.replace("38", "37")
        grch37_gff = str(Path(grch37_dir)/grch37_gff_name)
        if not os.path.isfile(grch37_gff):
            self.convert_38_to_37(grch38_gff, grch37_gff)
        # os.remove(grch38_gff)
        self._ftp_cnn.quit()

        return grch38_gff

    def convert_38_to_37(self, input_gff, output_gff) -> str:
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
