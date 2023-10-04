import os
import sys
import subprocess
import logging
logger = logging.getLogger(__name__)

class InstallationError(Exception):
    pass

def check_crossmap():
    '''
        Chheck if CrossMap.py is available on PATH
    '''
    cmd = "CrossMap.py -h"
    try:
        p1 = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    except:
        msg = (" ERROR: cannot execute command {}." \
            "Please, install CrossMap and make it accessible through the PATH").format(cmd)
        logging.error(msg)
        raise InstallationError(msg)
    else:
        output = p1.stdout.decode('UTF-8')
        error  = p1.stderr.decode('UTF-8')
        print(error)
        if error:
            msg = (" ERROR: cannot execute command {}." \
                "Please, install CrossMap and make it accessible through the PATH").format(cmd)
            logging.error(msg)
            raise InstallationError(error)
        else:
            msg = ("INFO: CrossMap.py is available on PATH")
            logging.info(msg)
