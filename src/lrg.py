import os
import sys
from pathlib import Path
from ftplib import FTP
class LRG():
    def __init__(self):
        pass
    @property
    def current_release(self):
        '''
        '''
        ftp = FTP("ftp.ebi.ac.uk")
        ftp.login()

        ftp.cwd("pub/databases/lrgex/")

        ori_file  = "current_schema_version.txt"
        dest_file = ori_file

        if not os.path.isfile(dest_file):
            ftp.retrbinary("RETR " + ori_file, open(dest_file, 'wb').write)

        latest_lrg_version = ""
        with open (dest_file) as f:
            for line in f:
                line = line.rstrip("\n")
                latest_lrg_version = line
        os.remove(dest_file)
        return latest_lrg_version

    def download_transcripts(self, output_dir):
        '''
            Download LRG transcripts with external references
        '''
        ftp = FTP("ftp.ebi.ac.uk")
        ftp.login()
        ftp.cwd("pub/databases/lrgex/")

        ori_file  = "list_LRGs_transcripts_xrefs.txt"
        dest_file = output_dir + "/" + ori_file

        if not os.path.isfile(dest_file):
            msg = ("INFO: Downloading {}").format(ori_file)
            print(msg)
            ftp.retrbinary("RETR " + ori_file, open(dest_file, 'wb').write)
        else:
            msg = ("INFO: File {} is already available at {}").format(ori_file, output_dir)
            print(msg)
        ftp.quit()
        return dest_file
