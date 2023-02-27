import argparse
from astropy.time import Time
from datetime import datetime
import json
import os
import requests
import pandas as pd
import subprocess
from lxml import html
import numpy as np
import re


# groupids  41:rcf, 43:clu
def api(method, endpoint, data=None):
    with open('/Users/viraj/ztf_utils/secrets_fritz_token.json', 'r') as f:
        dat = json.load(f)
    token = dat['token']
    headers = {'Authorization': f'token {token}'}
    response = requests.request(method, endpoint, params=data, headers=headers)
    return response.json()


def submit_fps_request(ra, dec, jd_start, jd_end):
    with open('/Users/viraj/ztf_utils/secrets_ztffps.json', 'r') as f:
        secrets = json.load(f)

    email = secrets['username']
    password = secrets['password']
    command = "wget --http-user=ztffps --http-passwd=dontgocrazy! -O log.txt 'https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?ra=%s&dec=%s&jdstart=%s&jdend=%s&email=%s&userpass=%s'" % (
    ra, dec, jd_start, jd_end, email, password)
    os.system(command)
    print(command)


def download_ztf_fphot(logdf, savefolder = 'fps_lcs'):
    '''
    Credit - Yashvi Sharma (lb)
    '''

    with open('/Users/viraj/ztf_utils/secrets_ztffps.json', 'r') as f:
        secrets = json.load(f)
    Email = secrets['username']
    UserPass = secrets['password']
    IUN = 'ztffps'
    IUW = 'dontgocrazy!'

    ldf = logdf[(logdf['service'] == 'ztf') & (logdf['status'] == 'pending')]
    if len(ldf) == 0:
        print('No pending requests, returning')
        return logdf
    inds = ldf.index
    statuspending = True
    StatusURL = 'https://ztfweb.ipac.caltech.edu/cgi-bin/getForcedPhotometryRequests.cgi'
    output = requests.get(StatusURL, auth=(IUN, IUW),
                          params={'email': Email, 'userpass': UserPass, 'option': 'All recent jobs',
                                  'action': 'Query Database'})  # past 30 days
    ohtml = html.fromstring(output.content)
    if (len(ohtml.xpath('//table')) == 0):
        print('No table, download manually')
        return logdf
    if (ohtml.xpath('//table')[0].text_content() == '') or (ohtml.xpath('//table')[0].text_content() is None):
        print('Empty table returned, no requests in last 30 days?')
        return logdf
    else:
        cols = [ele.text_content().strip() for ele in ohtml.xpath('//table/tr/th')]
        data = [re.split('[ ]{2,}', ele.text_content())[1:] for ele in ohtml.xpath('//table/tr')[1:]]
        outtbl = pd.DataFrame(columns=cols, data=data)
        for ind in inds:
            logrow = logdf.loc[ind]
            lcpath = list(outtbl[np.isclose(np.array(outtbl['ra']).astype(float), float(logrow['ra']),
                                            atol=1e-3) & np.isclose(np.array(outtbl['dec']).astype(float),
                                                                    float(logrow['dec']), atol=1e-3) & (
                                             outtbl['ended'] != '')]['lightcurve'])
            if len(lcpath) == 0:
                print('Job not finished for ', logrow['ztfname'])
                continue
            elif lcpath[0] is None:
                print(outtbl['ended'])
                print('LC path is None?')
                continue
            else:
                ztffphot_folder = savefolder
                url = f'https://ztfweb.ipac.caltech.edu{lcpath[0]}'
                fname = url.split('/')[-1]
                ztfname = logrow['ztfname']
                subprocess.call(
                    f'wget -P "{ztffphot_folder}"' + ' --http-user="ztffps" --http-passwd="dontgocrazy!" ' + url,
                    shell=True)
                # ztffphot_folder = dirpath + 'Work/ZTF/forcedphotometry_lightcurves/ztf_fphot'
                subprocess.call(f'mv {ztffphot_folder}/{fname} {ztffphot_folder}/{ztfname}_fps.dat',
                                shell=True)
                logdf.loc[ind, 'status'] = 'completed'
                print('Data downloaded for ', ztfname)
        return logdf


if __name__ == '__main__':
    today = Time(datetime.utcnow()).jd
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('action', choices=['request', 'plot', 'download'])
    parser.add_argument("name", type=str, help="Name of the target")
    parser.add_argument("--ra", type=str, help="RA of the target (decimal degrees)")
    parser.add_argument("--dec", type=str, help="Dec of the target (decimal degrees)")
    parser.add_argument("--fritz", action="store_true", help="Query Fritz for coordinates")
    parser.add_argument("--jd_start", type=float, default=2458000, help="start date in JD")
    parser.add_argument("--jd_end", type=float, default=today, help="end date in JD, default=today")
    parser.add_argument("--requests_logfile", type=str, default='request_log.csv', help="request logfile")

    args = parser.parse_args()
    request_logfile = args.requests_logfile

    if args.action == 'request':
        if args.fritz:
            response = api('GET', 'https://fritz.science/api/sources/%s' % (args.name))
            source = response['data']
            ra = source['ra']
            dec = source['dec']
            print(ra, dec)


        else:
            ra = args.ra
            dec = args.dec

        submit_fps_request(ra, dec, args.jd_start, args.jd_end)

        if not os.path.exists(request_logfile):
            logdf = pd.DataFrame({'ztfname': [args.name], 'ra': [ra],
                                  'dec': [dec], 'url': [None], 'status': ['pending'], 'service': ['ztf']})
            logdf.to_csv(request_logfile, index=False)

        else:
            logdf = pd.read_csv(request_logfile)
            newdf = pd.DataFrame({'ztfname': [args.name], 'ra': [ra],
                                  'dec': [dec], 'url': [None], 'status': ['pending'], 'service': ['ztf']})
            logdf = pd.concat([logdf,newdf])
            logdf.to_csv(request_logfile, index=False)

    if args.action == 'download':
        logdf = pd.read_csv(request_logfile)
        logdf = download_ztf_fphot(logdf)
        logdf.to_csv(request_logfile, index=False)