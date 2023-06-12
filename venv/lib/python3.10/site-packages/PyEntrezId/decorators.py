# coding: utf-8
import requests
import sys


hgnc = 9245
bad_hgnc = 'ENST'
good_server_call = "http://rest.genenames.org/fetch/hgnc_id/{0}".format(hgnc)
bad_server_call = "http://rest.genenames.org/jam/hgnc_id/{0}".format(bad_hgnc)
headers = {"Content-Type": "application/json"}


def retry(func):
    def retry_func(*args, **kwargs):
        max_tries = 3
        tries = 0
        while True:
            resp = func(*args, **kwargs)
            if not resp.ok and tries < max_tries:
                tries += 1
                print("Try # {0} - Please try again".format(tries))
                continue
            elif resp.ok:
                break
            else:
                resp.raise_for_status()
                sys.exit()
        return resp
    return retry_func


@retry
def get_response(server, headers):
    response = requests.get(server, headers=headers)
    return response
