from penquins import Kowalski
from astropy.io import ascii
import argparse
import json


def connect_kowalski():
    with open('secrets_kowalski.json','r') as f:
        secrets = json.load(f)

    username_kowalski = secrets['username']
    password_kowalski = secrets['password']
    # load username & passw credentials
    protocol, host, port = "https", "kowalski.caltech.edu", 443
    kowalski = Kowalski(username=username_kowalski, password=password_kowalski, protocol=protocol, host=host, port=port)
    connection_ok = kowalski.ping()
    print(f'Connection OK: {connection_ok}')
    return kowalski


def cone_search(k, coords, radius=20, catalog_name="ZTF_alerts"):
    # radius : arcsec
    # coords : {'ZTFname1:[ra1,dec1]'}
    q = {
        "query_type": "cone_search",
        "query": {
            "object_coordinates": {
                "cone_search_radius": radius,
                "cone_search_unit": "arcsec",
                "radec": coords
            },
            "catalogs": {
                catalog_name: {
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

    # print(q)
    r = k.query(query=q)
    data = r.get('data')
    return data


def cone_search_CLU(k, coords, radius=20, catalog_filter={}):
    # radius : arcsec
    # coords : {'ZTFname1:[ra1,dec1]'}
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
                    "filter": catalog_filter,
                    "projection": {
                        "_id": 0
                    }
                }
            }
        },
        "kwargs": {
            "filter_first": False
        }
    }

    # print(q)
    r = k.query(query=q)
    data = r.get('data')
    return data


def find_search(k, zmax):
    # Run a find query on CLU catalog on Kowalski for galaxies
    q = {
        'query_type': 'find',
        'query': {
            'catalog': 'CLU_20190625',
            'filter': {
                'z': {'$lt': zmax}
            },
            "projection": {
                "_id": 0,
                "name": 1,
                'z': 1,
                "ra": 1,
                "dec": 1,
                "magb": 1,
                "magberr": 1,
                "mstar": 1,
                "mstarerr": 1
            }
        }
    }

    r = k.query(query=q)
    data = r.get('data')
    return data


def find_search_aux(k, name):
    # Run a find query on CLU catalog on Kowalski for galaxies
    q = {
        'query_type': 'find',
        'query': {
            'catalog': 'ZTF_alerts_aux',
            'filter': {
                '_id': name
            },
            'projection': {
                'prv_candidates': 0
            },
        }
    }

    r = k.query(query=q)
    data = r.get('data')
    return data

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--name", type=str, help="object name")
    parser.add_argument("--ra", type=float, help="ra")
    parser.add_argument("--dec", type=float, help="dec")
    parser.add_argument("--radius", type=float, default=20, help="radius")
    parser.add_argument("--clu", action="store_true")
    parser.add_argument("--aux", action="store_true")
    parser.add_argument("--full", action="store_true",help="Print full ZTF alert")

    args = parser.parse_args()

    k = connect_kowalski()
    coords = {args.name: [args.ra, args.dec]}
    data = cone_search(k, coords, args.radius)

    if len(data['ZTF_alerts'][args.name]) > 0:
        for i in data['ZTF_alerts'][args.name]:
            print(i['objectId'])
            if args.full:
                print(i['candidate'])

    else:
        print(data)

    if args.clu:
        print('CLU')
        data = cone_search_CLU(k, coords, args.radius)
        if len(data['CLU_20190625'][args.name]) > 0:
            print(data['CLU_20190625'][args.name])
        else:
            print(data)

    if args.aux:
        print('ZTF alerts aux')
        # data = cone_search(k, coords=coords, radius=args.radius, catalog_name="ZTF_alerts_aux")
        data = find_search_aux(k, name=args.name)
        print(data)