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

    def __init__(self, genome_version):
        self._genome_version = genome_version

    @property
    def current_release(self)->str:
        '''
        '''
        ftp = FTP("ftp.ncbi.nlm.nih.gov")
        ftp.login()

        latest_folder = "refseq/MANE/MANE_human/current/"
        ftp.cwd(latest_folder)

        ori_file  = "README_versions.txt"
        dest_file = ori_file

        if not os.path.isfile(dest_file):
            ftp.retrbinary("RETR " + ori_file, open(dest_file, 'wb').write)

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

    def download_gff(self, output_dir) -> str:
        '''
        '''
        ftp = FTP("ftp.ncbi.nlm.nih.gov")
        ftp.login()

        latest_folder = "refseq/MANE/MANE_human/current/"
        ftp.cwd(latest_folder)

        grch38_dir =  output_dir + "/" + self.current_release + "/" + "GRCh38"
        Path(grch38_dir).mkdir(parents=True, exist_ok=True)

        grch38_gff_name  = "MANE.GRCh38.v{}.ensembl_genomic.gff.gz".format(self.current_release)
        grch38_gff = str(Path(grch38_dir)/grch38_gff_name)
        if not os.path.isfile(grch38_gff):
            ftp.retrbinary("RETR " + grch38_gff_name, open(grch38_gff, 'wb').write)

        grch37_dir =  output_dir + "/" + self.current_release + "/" + "GRCh37"
        Path(grch37_dir).mkdir(parents=True, exist_ok=True)
        grch37_gff_name = grch38_gff_name.replace("38", "37")
        grch37_gff = str(Path(grch37_dir)/grch37_gff_name)
        if not os.path.isfile(grch37_gff):
            self.convert_38_to_37(grch38_gff, grch37_gff)
        # os.remove(grch38_gff)

        return grch38_gff

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
