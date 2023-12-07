import numpy as np
import requests
import json
from astropy.time import Time
from datetime import datetime
import argparse
import pandas as pd
import os


def submit_batch_fps(ra_list, dec_list, jd_start=2458194.6,
                     jd_end=Time(datetime.utcnow()).jd):
    if not os.path.exists('secrets_ztffps.json'):
        raise FileNotFoundError("Please save your secrets in the file secrets_ztffps.json with the format {'username':<>, 'password': <>}")
    with open('secrets_ztffps.json', 'r') as f:
        secrets = json.load(f)['default_credentials']
    Email = secrets['username']
    UserPass = secrets['password']

    ra = json.dumps(ra_list)
    dec = json.dumps(dec_list)

    jdstart = json.dumps(jd_start)

    jdend = json.dumps(jd_end)
    print(jdend)
    payload = {'ra': ra, 'dec': dec,
               'jdstart': jdstart, 'jdend': jdend,
               'email': Email, 'userpass': UserPass}

    url = 'https://ztfweb.ipac.caltech.edu/cgi-bin/batchfp.py/submit'
    r = requests.post(url, auth=('ztffps', 'dontgocrazy!'), data=payload)
    print(r.status_code)


def retrieve_batch_fps(logfile='ztf_batch_fps_log.csv'):
    with open('secrets_ztffps.json', 'r') as f:
        secrets = json.load(f)['default_credentials']
    Email = secrets['username']
    UserPass = secrets['password']

    settings = {'email': Email,'userpass': UserPass,
                'option': 'All recent jobs', 'action': 'Query Database'}
    r = requests.get('https://ztfweb.ipac.caltech.edu/cgi-bin/getBatchForcedPhotometryRequests.cgi',
    auth=('ztffps', 'dontgocrazy!'),params=settings)
    df = pd.read_html(r.text)[0]
    df.to_csv(logfile)
    print(r.status_code)


def download_lightcurves(df, lcdir='batch_fps_lcs'):
    for ind in range(len(df)):
        row = df.loc[ind]
        reqid = row['reqId']
        lcname = f"{lcdir}/batchfp_{reqid}_lc.txt"

        depot_lcpath = row['lightcurve']
        os.system(f'''wget --http-user=ztffps --http-passwd=dontgocrazy! -O {lcname} "https://ztfweb.ipac.caltech.edu{depot_lcpath}"''')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['submit', 'retrieve', 'download'], help='Submit or retrieve?')
    parser.add_argument('-filename', type=str, help='Name of file to perform bfp on', default=None)
    parser.add_argument('-jd_start', default=2458194.6, help='JD Start')
    parser.add_argument('-jd_end', default=Time(datetime.utcnow()).jd, help='JD End')
    parser.add_argument('-ra_key', default='ra [deg]', help='RA key')
    parser.add_argument('-dec_key', default='dec [deg]', help='Dec key')
    args = parser.parse_args()

    if args.mode == 'submit':
        if args.filename is None:
            raise ValueError("Please provide -filename")

        bfp_filename = args.filename

        bfp_df = pd.read_csv(bfp_filename)

        bfp_df_split = np.array_split(bfp_df, int(len(bfp_df)/1500)+1)

        for subdf in bfp_df_split:
            submit_batch_fps(ra_list=list(subdf[args.ra_key]),
                             dec_list=list(subdf[args.dec_key]),
                             jd_start=args.jd_start,
                             jd_end=args.jd_end)

    if args.mode == 'retrieve':
        retrieve_batch_fps()

    if args.mode == 'download':
        if args.filename is None:
            raise ValueError("Please provide -filename")

        df = pd.read_csv(args.filename)

        download_lightcurves(df)
        

