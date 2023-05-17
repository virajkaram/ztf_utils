import os
import pandas as pd
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-names', type=str, nargs='+', default=None)
    parser.add_argument('-names_file', type=str, default=None)
    parser.add_argument('-time_isot_vline', type=str, default=None)
    args = parser.parse_args()

    if args.names is not None:
        names = args.names
        lcpaths = [f'fps_lcs/{name}_fps.dat' for name in names]
    if args.names_file is not None:
        names_df = pd.read_csv(args.names_file)
        names = names_df['ZTF_name']
        if 'lcpath' in names_df:
            lcpaths = names_df['lcpath']
        else:
            lcpaths = [f'fps_lcs/{name}_fps.dat' for name in names]

    for lcpath in lcpaths:
        cmd = f'python plot_fps.py --ztf_file {lcpath} ' \
              f'--mag --ymin 20 --ymax 23 --xmin -10'

        if args.time_isot_vline is not None:
            cmd += f' --time_isot_vline {args.time_isot_vline}'

        os.system(cmd)
