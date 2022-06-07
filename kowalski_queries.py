from penquins import Kowalski
from astropy.io import ascii
import argparse


def connect_kowalski():
    secrets = ascii.read('/Users/viraj/ztf_utils/secrets.csv', format = 'csv')
    username_kowalski = secrets['kowalski_user'][0]
    password_kowalski = secrets['kowalski_pwd'][0]
    # load username & passw credentials
    protocol, host, port = "https", "kowalski.caltech.edu", 443
    kowalski = Kowalski(username=username_kowalski, password=password_kowalski,protocol=protocol,host=host,port=port)
    connection_ok = kowalski.ping()
    print(f'Connection OK: {connection_ok}')
    return kowalski


def cone_search(k,coords, radius=20):
    #radius : arcsec
    #coords : {'ZTFname1:[ra1,dec1]'}
    q = {
        "query_type": "cone_search",
        "query": {
            "object_coordinates": {
                "cone_search_radius": radius,
                "cone_search_unit": "arcsec",
                "radec": coords
            },
            "catalogs": {
                "ZTF_alerts": {
                    "filter": {},
                    "projection": {
                        "_id": 0,
                        
                    }
                }
            }
        },
        "kwargs": {
            "filter_first": False
        }
    }
    
    #print(q)
    r = k.query(query=q)
    data = r.get('data')
    return data


def cone_search_CLU(k, coords, radius=20):
    #radius : arcsec
    #coords : {'ZTFname1:[ra1,dec1]'}
    q = {
        "query_type": "cone_search",
        "query": {
            "object_coordinates": {
                "cone_search_radius": radius,
                "cone_search_unit": "arcsec",
                "radec": coords
            },
            "catalogs": {
                "CLU_20190625": {
                    "filter": {},
                    "projection": {
                        "_id": 0,
                        
                    }
                }
            }
        },
        "kwargs": {
            "filter_first": False
        }
    }
    
    #print(q)
    r = k.query(query=q)
    data = r.get('data')
    return data


def find_search(zmax):
    #Run a find query on CLU catalog on Kowalski for galaxies
    q = {
    'query_type': 'find',
    'query': {
        'catalog': 'CLU_20190625',
        'filter': {
            'z': {'$lt':zmax}
        },
        "projection": {
            "_id": 0,
            "name": 1,
            'z': 1,
            "ra":1,
            "dec":1,
            "magb":1,
            "magberr":1,
            "mstar":1,
            "mstarerr":1
        }
    }
    }


    r = k.query(query=q)
    data = r.get('data')
    return data


if __name__ =='__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--name",type=str,help="object name")
    parser.add_argument("--ra",type=float,help="ra")
    parser.add_argument("--dec",type=float,help="dec")
    parser.add_argument("--radius",type=float,default=20,help="radius")
    parser.add_argument("--clu",action="store_true")

    args = parser.parse_args()

    k = connect_kowalski()
    coords = {args.name:[args.ra,args.dec]}
    data = cone_search(k,coords,args.radius)
    if len(data['ZTF_alerts'][args.name])>0:
        for i in data['ZTF_alerts'][args.name]:
            print(i['objectId'])

    else:
        print(data)

    if args.clu:
        print('CLU')
        data = cone_search_CLU(k,coords,args.radius)
        if len(data['CLU_20190625'][args.name]) > 0:
            print(data['CLU_20190625'][args.name])
        else:
            print(data)
