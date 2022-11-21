from utils import in_ellipse, api
import argparse
from kowalski_queries import cone_search_CLU, connect_kowalski
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.units as u


def get_xmatched_galaxies(k, ra, dec, name, radius_kpc=30):
    radius_deg = 3
    radius_arcsec = 3 * 3600
    response = cone_search_CLU(k=k, coords={name: [ra, dec]}, radius=radius_arcsec)
    galaxies = response['CLU_20190625'][name]
    print(f'Found {len(galaxies)} galaxies in coarse search.')
    # these guys are very big, so check them separately
    M31 = {
        "_id": 596900,
        "name": "PGC2557",
        "ra": 10.6847,
        "dec": 41.26901,
        "a": 6.35156,
        "b2a": 0.32,
        "pa": 35.0,
        "z": -0.00100100006,
        "sfr_fuv": None,
        "mstar": 253816876.412914,
        "sfr_ha": 0,
        "coordinates": {"radec_str": ["00:42:44.3503", "41:16:08.634"]},
    }
    M33 = {
        "_id": 597543,
        "name": "PGC5818",
        "ra": 23.46204,
        "dec": 30.66022,
        "a": 2.35983,
        "b2a": 0.59,
        "pa": 23.0,
        "z": -0.000597000006,
        "sfr_fuv": None,
        "mstar": 4502777.420493,
        "sfr_ha": 0,
        "coordinates": {"radec_str": ["01:33:50.8900", "30:39:36.800"]},
    }

    matches = []
    for galaxy in galaxies + [M31, M33]:
        alpha1, delta01 = galaxy["ra"], galaxy["dec"]
        redshift = galaxy["z"]

        if redshift < 0.01:
            # for nearby galaxies and galaxies with negative redshifts, do a 5 arc-minute cross-match
            # (cross-match radius would otherwise get un-physically large for nearby galaxies)
            cm_radius = 300.0 / 3600

        else:
            cm_radius = radius_kpc * (0.05 / redshift) / 3600

        in_galaxy = in_ellipse(ra, dec, alpha1, delta01, cm_radius, 1, 0)
        if in_galaxy:
            match = galaxy
            matches.append(match)

    return matches


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--ra", type=float, help="object ra")
    parser.add_argument("--dec", type=float, help="object dec")
    parser.add_argument("--name", type=str, help="Name", default='ZTF')
    parser.add_argument("--fritz", action="store_true", help="Query Fritz for coordinates")
    parser.add_argument("--file", type=str, help='Names in a file')
    # parser.add_argument("--gra", type=float, help="galaxy ra")
    # parser.add_argument("--gdec", type=float, help="galaxy dec")
    # parser.add_argument("--gz", type=float, help="galaxy redshift")

    args = parser.parse_args()
    # alpha1, delta01, redshift = args.gra, args.gdec, args.gz

    k = connect_kowalski()

    if args.file is not None:
        transients = ascii.read(args.file)
        crds = SkyCoord(ra=transients['RA'], dec=transients['Dec'], unit=(u.hour, u.deg))

        matched_names = []
        for ind in range(len(transients)):
            crd = crds[ind]
            name = transients[ind]['ZTFID']
            ra, dec = crd.ra.deg, crd.dec.deg
            matches = get_xmatched_galaxies(k, ra, dec, name=name)
            galaxy_names = [x['name'] for x in matches]
            matched_names.append(galaxy_names)
            print(f"{ind}/{len(transients)}, {galaxy_names}")

        with open('local_rcf_not_in_clu_clu_xmatch.tsv', 'w') as f:
            f.write(f'Name\tClu_match\n')
            for ind in range(len(transients)):
                f.write(f"{transients[ind]['ZTFID']}\t{matched_names[ind]}\n")

    else:
        if args.fritz:
            name = args.name
            response = api('GET', f'https://fritz.science/api/sources/{name}')
            source = response['data']
            ra = source['ra']
            dec = source['dec']
            print(ra, dec)

        else:
            ra, dec, name = args.ra, args.dec, args.name

        matches = get_xmatched_galaxies(k, ra, dec, name=name)
        names = [x['name'] for x in matches]
        print(f"Cross-matched to {names}")
