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


def submit_fps_request(ra,dec,jd_start,jd_end):
	with open('/Users/viraj/ztf_utils/secrets_ztffps.json','r') as f:
		secrets = json.load(f)

	email = secrets['username']
	password = secrets['password']
	command = "wget --http-user=ztffps --http-passwd=dontgocrazy! -O log.txt 'https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?ra=%s&dec=%s&jdstart=%s&jdend=%s&email=%s&userpass=%s'"%(ra,dec,jd_start,jd_end,email,password)
	os.system(command)
	print(command)



if __name__ == '__main__':
    today = Time(datetime.utcnow()).jd
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("name",type=str,help="Name of the target")
    parser.add_argument("--ra",type=str,help="RA of the target (decimal degrees)")
    parser.add_argument("--dec",type=str,help="Dec of the target (decimal degrees)")
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


    submit_fps_request(ra,dec,args.jd_start,args.jd_end)
