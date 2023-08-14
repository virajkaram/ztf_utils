import numpy as np
import pandas as pd
import os
from glob import glob
import requests, json
from astropy.io import ascii, fits
from astropy.time import Time
from astropy.table import Table
import matplotlib.pyplot as plt
import argparse


def api(method, endpoint, data=None):
    ''' Info : Basic API query, takes input the method (eg. GET, POST, etc.), the endpoint (i.e. API url)
               and additional data for filtering
        Returns : response in json format
        CAUTION! : If the query doesn't go through, try putting the 'data' input in 'data' or 'params'
                    argument in requests.request call
    '''
    headers = {'Authorization': f'token {token}'}
    response = requests.request(method, endpoint, json=data, headers=headers)
    return response.json()


def get_instrument_id(inst_name):
    url = BASEURL + 'api/instrument'
    response = api('GET', url)
    inst_list = np.array(response['data'])
    inst_id = inst_list[[inst_name in d['name'] for d in inst_list]][0]['id']
    return inst_id


def get_user_id(firstname, lastname, loosematch=False):
    url = BASEURL + 'api/user'
    response = api('GET', url)
    user_list = response['data']['users']
    user_info = []
    for d in user_list:
        fn = 'NA' if (d['first_name'] is None) else d['first_name']
        ln = 'NA' if (d['last_name'] is None) else d['last_name']
        if (loosematch):
            if (firstname in fn or lastname in ln):
                user_info.append(d)
        else:
            if (firstname in fn and lastname in ln):
                user_info.append(d)
                break
    if (len(user_info) > 1):
        print('Following users found, choose array index of your user')
        print(user_info)
        inp = input('Array index? : ')
        user = user_info[int(inp)]
    else:
        print('User found', user_info)
        user = user_info[0]
    return user['id']


def get_source_api(ztfname):
    ''' Info : Query a single source, takes input ZTF name
        Returns : all basic data of that source (excludes photometry and spectra,
                  includes redshift, classification, comments, etc.)
    '''
    url = BASEURL + 'api/sources/' + ztfname
    response = api('GET', url)
    return response['data']


def upload_spectrum_from_ascii(ztfname, instrument_id, filename, obsdate, filedata, observers, reducers, groups):
    '''Info : Uploads ascii filedata to ztfname, NEEDS instrument_id, obsdate(UTC) and filename
       Returns : Response with status
    '''
    url = BASEURL + 'api/spectrum/ascii'
    if instrument_id == 7:
        filt = {'obj_id': ztfname, 'instrument_id': instrument_id, 'filename': filename, 'observed_at': str(obsdate),
                'ascii': filedata, 'observed_by': observers, 'reduced_by': reducers, 'group_ids': groups,
                'wave_column': 0, 'flux_column': 1, 'fluxerr_column': 3}
    elif instrument_id == 3:
        filt = {'obj_id': ztfname, 'instrument_id': instrument_id, 'filename': filename, 'observed_at': str(obsdate),
                'ascii': filedata, 'observed_by': observers, 'reduced_by': reducers, 'group_ids': groups,
                'wave_column': 0, 'flux_column': 1, 'fluxerr_column': 2}
    response = api('POST', url, filt)
    return response


def upload_spec(specfile, ztfname, obsdate, meta, instrument_id, observers, reducers):
    filedata = open(specfile, 'r').read()
    print(ztfname, get_source_api(ztfname))
    sourcegroups = get_source_api(ztfname)['groups']
    groupids = [sg['id'] for sg in sourcegroups]
    spec = specfile
    print(ztfname, obsdate)
    inp = input('upload?')
    if inp == 'y':
        r = upload_spectrum_from_ascii(ztfname, instrument_id, os.path.basename(specfile), obsdate, filedata, observers,
                                       reducers, groupids)
        print(r)


def prepare_spectrum_in_ascii_format(spec, instrument_id, clip=True, minwav=3500):
    ascii_specname = ''
    if instrument_id == 3:
        name = os.path.basename(spec).split('_')[0]
        hdr = fits.open(spec)[1].header
        dat = fits.open(spec)[-1].data
        tbl = Table(dat)
        tbl = tbl[tbl['wave'] > 3500]
        hdr = dict(hdr)
        hdr.pop('COMMENT')
        tbl.meta['comments'] = [key + ': ' + str(hdr[key]) for key in list(hdr.keys())] + [
            '# Columns: WAVE FLUX FLUX_ERR']
        ascii_specname = f'{args.d}/{name}_spec.txt'
        ascii.write(tbl, ascii_specname, format='no_header', overwrite=True)
        lines = np.array(tbl.meta['comments'])
        date = lines[['UTSHUT' in line for line in lines]][0].split(': ')[1]

    if instrument_id==7:
        lines = open(spec, 'r').readlines()
        date = Time(np.array(lines)[['MJD' in line for line in lines]][0].split()[3][1:-1], format='mjd').isot
        if clip:
            tbl = ascii.read(spec)
            tbl = tbl[tbl['col1']>minwav]
            ascii.write(tbl, spec, format='no_header', overwrite=True)
        ascii_specname = spec

    return ascii_specname, date


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str, help='Directory')
    parser.add_argument('-obs_names', type=str,
                        help='Comma separated list of observer names e.g Viraj Karambelkar, Nicholas Earley')
    parser.add_argument('-red_names', type=str, help='Comma separated list of reducer names e.g. Viraj Karambelkar')
    parser.add_argument('-inst', type=str, default='DBSP', help='Instrument name, e.g. DBSP, LRIS')
    args = parser.parse_args()

    global token, BASEURL
    token = os.getenv('FRITZ_TOKEN')  ## generate token from your profile and add here
    if token is None:
        print(f"Fritz token not set. Please run 'export FRITZ_TOKEN=<your_fritz_token>' ")
        exit(1)

    BASEURL = 'https://fritz.science/'

    instrument_id = get_instrument_id(args.inst)
    observers = args.obs_names.split(',')
    observer_ids = []
    for observer in observers:
        firstname, lastname = observer.split(' ')
        observer_ids.append(get_user_id(firstname, lastname))

    print(f'Observer ids : {observer_ids}')
    reducers = args.red_names.split(',')
    reducer_ids = []
    for reducer in reducers:
        firstname, lastname = reducer.split(' ')
        reducer_ids.append(get_user_id(firstname, lastname))
    print(f'Reducer ids : {reducer_ids}')
    
    if instrument_id==3:    
        speclist = glob(f'{args.d}/*.fits')
    if instrument_id == 7:
        speclist = glob(f'{args.d}/*.spec')
    ascii_speclist = []
    uploaded = []

    for spec in speclist:
        if instrument_id==3:
            ztfname = os.path.basename(spec).split('_')[0]
        if instrument_id==7:
            ztfname = os.path.basename(spec).split('_')[1]
        print(ztfname)
        if ztfname in uploaded:
            continue
        ascii_specname, date = prepare_spectrum_in_ascii_format(spec, instrument_id)
        spec = ascii.read(ascii_specname)
        plt.figure()
        plt.plot(spec['col1'], spec['col2'])
        plt.show()
        upload_spec(specfile=ascii_specname, ztfname=ztfname, obsdate=date, meta=None, instrument_id=instrument_id,
                    observers=observer_ids, reducers=reducer_ids)
        uploaded.append(ztfname)

    print(f"Successfully uploaded {uploaded} : {len(uploaded)}/{len(speclist)}")
