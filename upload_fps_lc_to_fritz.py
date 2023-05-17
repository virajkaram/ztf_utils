import base64
import pandas as pd
import argparse
import requests
import json


def pdf_to_base64(filepath):
    with open(filepath, "rb") as pdf_file:
        encoded_string = base64.b64encode(pdf_file.read())
    return encoded_string


def api(method, endpoint, data=None):
    with open('/Users/viraj/ztf_utils/secrets_fritz_token.json', 'r') as f:
        dat = json.load(f)
    token = dat['token']
    headers = {'Authorization': f'token {token}'}
    print(method, endpoint, headers, data)
    response = requests.request(method, endpoint, params=data, headers=headers)
    return response.json()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-names', type=str, nargs='+', default=None)
    parser.add_argument('-names_file', type=str, default=None)

    args = parser.parse_args()

    if args.names is not None:
        names = args.names
        lcpaths = [f'fps_lcs/{name}_fps_mags.pdf' for name in names]

    if args.names_file is not None:
        names_df = pd.read_csv(args.names_file)
        names = names_df['ZTF_name']
        if 'lcpath' in names_df:
            lcpaths = names_df['lcpath']
        else:
            lcpaths = [f'fps_lcs/{name}_fps_mags.pdf' for name in names]

    for name, lcpath in zip(names, lcpaths):

        base64_lc = pdf_to_base64(lcpath)

        data = {'text': 'test comment',
                # 'attachment': {'body': base64_lc, 'name': f'{name}_fps.dat'}
                }

        print(f"Uploading {name}")
        response = api('POST', f'https://fritz.science/api/sources/{name}/comments',
                       data=data)
        print(response)
