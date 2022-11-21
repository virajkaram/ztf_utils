import argparse
from penquins import Kowalski
from astropy.io import ascii
import json
from kowalski_queries import connect_kowalski

def query_aux_Kowalski(k, name):
    #Query CLU catalog on Kowalski for galaxies around coords within a distance of xx arcsec
    #coords = {'ZTFname1':[ra1,dec1],'ZTFname2':[ra2,dec2]}
    q = {
    'query_type': 'find',
    'query': {
        'catalog': 'ZTF_alerts_aux',
        'filter': {
            '_id': name
        },
        "projection": {
            "_id": 0
        }
    }
    }

    r = k.query(query=q)
    data = r.get('data')
    return data

def query_alerts_Kowalski(k, name):
    #Query CLU catalog on Kowalski for galaxies around coords within a distance of xx arcsec
    #coords = {'ZTFname1':[ra1,dec1],'ZTFname2':[ra2,dec2]}
    q = {
    'query_type': 'find',
    'query': {
        'catalog': 'ZTF_alerts',
        'filter': {
            'objectId': name
        },
        "projection": {
            "_id": 0,
            "cutoutScience":0,
            "cutoutTemplate":0,
            "cutoutDifference":0
        }
    }
    }

    r = k.query(query=q)
    data = r.get('data')
    return data


def save_alerts(k, name):
    aux_data = query_aux_Kowalski(k, name=name)
    alerts = query_alerts_Kowalski(k, name=name)
    for ind, alert in enumerate(alerts):
        alert['cross_matches'] = aux_data[0]['cross_matches']
        alert['prv_candidates'] = aux_data[0]['prv_candidates']
        with open(f'alerts/alert_{name}_{ind}.json', 'w') as f:
            json.dump([alert], f)


if __name__ =='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--name',type=str, help='ZTF name')

    args = parser.parse_args()
    name = args.name

    k = connect_kowalski()
    save_alerts(k,name)