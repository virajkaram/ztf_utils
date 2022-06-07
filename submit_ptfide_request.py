import argparse
from astropy.time import Time
from datetime import datetime
import json
import os
import requests

#groupids  41:rcf, 43:clu
def api(method, endpoint, data=None):
    with open('/Users/viraj/ztf_utils/secrets_fritz_token.json','r') as f:
        dat = json.load(f)
    token = dat['token']
    headers = {'Authorization': f'token {token}'}
    response = requests.request(method, endpoint, params=data, headers=headers)
    return response.json()


def submit_ptfide_request(ra,dec,ptffieldid,fid,ccdid,jd_start,jd_end,name):
    with open('/Users/viraj/ztf_utils/secrets_ptffps.json','r') as f:
        dat = json.load(f)
    user = dat['username']
    pwd = dat['password']
    command = "wget --http-user=%s --http-passwd=%s -O query.html 'http://ptfdepot.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?ra=%s&dec=%s&ptffield=%s&fid=%s&ccdid=%s&jdstartoflightcurve=%s&jdendoflightcurve=%s&email=viraj.karambelkar@gmail.com&priority=1&name=%s'"%(user,pwd,ra,dec,ptffieldid,fid,ccdid,jd_start,jd_end,name)
    print(command)
    os.system(command)



if __name__ == '__main__':
    today = Time(datetime.utcnow()).jd
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("name",type=str,help="Name of the target")
    parser.add_argument("--ra",type=str,help="RA of the target (decimal degrees)")
    parser.add_argument("--dec",type=str,help="Dec of the target (decimal degrees)")
    parser.add_argument("--fid",type=int,help="Filter ID")
    parser.add_argument("--field",type=int,help="Field ID")
    parser.add_argument("--ccdid",type=int,help="CCD ID")
    parser.add_argument("--fritz",action="store_true",help="Query Fritz for coordinates")
    parser.add_argument("--jd_start",type=float,default=2458000,help="start date in JD")
    parser.add_argument("--jd_end",type=float,default=today,help="end date in JD, default=today")

    args = parser.parse_args()
    token = 'd192b631-d073-4c1a-85d0-3a8e14a130f7'

    if args.fritz:
    	response = api('GET', 'https://fritz.science/api/sources/%s'%(args.name))
    	source = response['data']
    	ra = source['ra']
    	dec = source['dec']
    	print(ra,dec)


    else :
    	ra = args.ra
    	dec = args.dec


    submit_ptfide_request(ra,dec,ptffieldid=args.field, fid=args.fid, ccdid=args.ccdid,jd_start=args.jd_start,jd_end=args.jd_end,name=args.name)
