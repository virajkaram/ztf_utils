import os
import pandas as pd
import argparse
from astropy.coordinates import SkyCoord
import astropy.units as u
import shutil

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str,
                        help='Name of file to perform bfp on', default=None)
    parser.add_argument('-lcdir', type=str,
                        help='Name of file to perform bfp on', default='lcs')
    args = parser.parse_args()

    df = pd.read_csv(args.filename)

    clu_sources = pd.read_csv('all_clu_sources_20230412_fps.csv')
    clu_sources_crds = SkyCoord(ra=clu_sources['ra [deg]'],
                                dec=clu_sources['dec [deg]'],
                                unit=(u.deg, u.deg))
    for ind, row in df.iterrows():
        ra, dec = row['ra'], row['dec']
        reqid = row['reqId']
        lcfile = f"{args.lcdir}/batchfp_{reqid}_lc.txt"

        crd = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
        idx, d2d, d3d = crd.match_to_catalog_sky(clu_sources_crds)
        ztf_name = clu_sources.iloc[idx]['id']
        ztf_name_lc = f"named_lcs/{ztf_name}_fps.dat"
        shutil.copy(lcfile, ztf_name_lc)
        os.system(f"python plot_fps.py --ztf_file {ztf_name_lc} --no_mask --mag")
        os.system(f"python plot_fps.py --ztf_file {ztf_name_lc} --no_mask --bin "
                  f"--rbin 3 --gbin 3 --ibin 3")
        os.system(f"")
