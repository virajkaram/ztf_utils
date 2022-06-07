import os
from astropy.io import ascii
from astropy.table import Table, Column 
import requests
import json

sourcenames = ascii.read('lrn-ilrt-ZTFI.csv')

def api(method, endpoint, data=None):
    with open('/Users/viraj/ztf_utils/secrets_fritz_token.json','r') as f:
        dat = json.load(f)
    token = dat['token']
    headers = {'Authorization': f'token {token}'}
    response = requests.request(method, endpoint, params=data, headers=headers)
    return response.json()

ras = []
decs = []
for source in sourcenames:
        print('Currently doing %s'%(source['Name']))
        response = api('GET', 'https://fritz.science/api/sources/%s'%(source['Name']))
        source = response['data']
        ra = source['ra']
        dec = source['dec']
        ras.append(ra)
        decs.append(dec)

sourcenames.add_column(Column(name='RA',data=ras))
sourcenames.add_column(Column(name='Dec',data=decs))

sourcenames.write('lrn-ilrt-ZTFI.csv',overwrite=True)
