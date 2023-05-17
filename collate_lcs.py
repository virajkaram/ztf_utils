import pandas as pd
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-names', type=str, nargs='+', default=None)
    parser.add_argument('-names_file', type=str, default=None)

    args = parser.parse_args()

    if args.names is not None:
        names = args.names
        lcpaths = [f'fps_lcs/{name}_fps_mags.pdf' for name in names]
        lcs_dir = 'fps_lcs_collection'

    if args.names_file is not None:
        names_df = pd.read_csv(args.names_file)
        names = names_df['ZTF_name']
        if 'lcpath' in names_df:
            lcpaths = names_df['lcpath']
        else:
            lcpaths = [f'fps_lcs/{name}_fps_mags.pdf' for name in names]

        lcs_dir = args.names_file.replace('.txt', '') + '_fps_lcs_collection'

    if os.path.exists(lcs_dir):
        os.system(f'rm -rf {lcs_dir}')
    os.mkdir(lcs_dir)
    for lcpath in lcpaths:
        os.system(f'cp {lcpath} {lcs_dir}/')
