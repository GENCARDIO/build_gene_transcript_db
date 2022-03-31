import os
import sys
from pathlib import Path
from ftplib import FTP
import logging
logger = logging.getLogger(__name__)

class LRG():
    '''
        Class that handles LRG transcripts human GFF data
        :param str output_dir: output directory
        :py:attr:`current_version`
        :py:meth:`download_transcripts`
    '''

    def __init__(self, output_dir):
        '''
        '''

        self._output_dir     = output_dir
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

        self._working_dir    = os.path.join(self._output_dir, "LRG",
            self.current_release)
        Path(self._working_dir).mkdir(parents=True, exist_ok=True)

    @property
    def current_release(self):
        '''
        '''
        self._ftp_cnn.cwd(self._parent_url)

        url = "pub/databases/lrgex/"
        try:
            self._ftp_cnn.cwd(url)
        except:
            msg = ("ERROR: Unable to move to {}").format(url)
            logging.error(msg)
            raise ConnectionError(msg)

        ori_file  = "current_schema_version.txt"
        dest_file = ori_file

        if not os.path.isfile(dest_file):
            self._ftp_cnn.retrbinary("RETR " + ori_file, open(dest_file, 'wb').write)

        latest_lrg_version = ""
        with open (dest_file) as f:
            for line in f:
                line = line.rstrip("\n")
                latest_lrg_version = line
        f.close()
        os.remove(dest_file)
        return latest_lrg_version

    def download_transcripts(self):
        '''
            Download LRG transcripts with external references
        '''
        self._ftp_cnn.cwd(self._parent_url)

        url = "pub/databases/lrgex/"
        try:
            self._ftp_cnn.cwd(url)
        except:
            msg = ("ERROR: Unable to move to {}").format(url)
            logging.error(msg)
            raise ConnectionError(msg)

        ori_file  = "list_LRGs_transcripts_xrefs.txt"
        dest_file = os.path.join(self._working_dir, ori_file)

        if not os.path.isfile(dest_file):
            msg = ("INFO: Downloading {}").format(ori_file)
            logging.error(msg)
            self._ftp_cnn.retrbinary("RETR " + ori_file, open(dest_file, 'wb').write)
        else:
            msg = ("INFO: File {} is already available at {}").format(ori_file, self._working_dir)
            logging.error(msg)
        self._ftp_cnn.quit()
        return dest_file
