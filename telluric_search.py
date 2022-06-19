import lxml.html
import requests
from astropy.io import ascii
import numpy as np
import argparse

def query_gemini_website(ra_hms, dec_dms, min_dist_deg=15, name='Science star', min_mag=8, max_mag=11):
    webpage = requests.get('https://www.gemini.edu/apps/telluric-search/telluric-search.py')
    doc = lxml.html.fromstring(webpage.content)
    form = doc.xpath('//form')[0]

    form.fields['ra'] = ra_hms
    form.fields['dec'] = dec_dms
    form.fields['sptype'] = 'A0V'
    form.action = 'https://www.gemini.edu/apps/telluric-search/telluric-search.py'

    olink = lxml.html.submit_form(form)
    opage = lxml.html.parse(olink).getroot()
    ochildren = opage.getchildren()
    opre = ochildren[1].getchildren()[0].find('pre')
    table_content = opre.text_content()
    table_content.replace('Simbad', '')
    telluric_stars = ascii.read(table_content.replace('Simbad', ''), delimiter=' ')

    mask = (telluric_stars['Dist(deg)'] < min_dist_deg) & (telluric_stars['Vmag'] < max_mag) & (
                telluric_stars['Vmag'] > min_mag)

    if np.all(telluric_stars['Vmag']<min_mag):
        print(
            f'All telluric stars for {name} are brighter than {min_mag} magnitudes. Returning the faintest among those')
        return telluric_stars[np.argmax(telluric_stars['Vmag'])]

    if len(telluric_stars[mask])==0:
        print(f'Not found any stars within the set distance and magnitude range. Trying with a larger magnitude range')
        mask = (telluric_stars['Dist(deg)'] < 15) & (telluric_stars['Vmag'] < max_mag) & (
                telluric_stars['Vmag'] > 5)
        if len(telluric_stars[mask])==0:
            print(f'Failed for {name}, no stars found')
            return -99
        else:
            print(f'Found a star in increased parameter search, please check manually to make sure it is okay')
    return telluric_stars[mask][0]


if __name__ =='__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input_filename", type=str, help='Input targetlist in Keck/P200 format')

    args = parser.parse_args()
    inp_filename = args.input_filename

    with open(inp_filename, 'r') as f:
        dat = f.readlines()

    sci_targets = []
    for row in dat:
        if np.logical_or('_o' in row, 'HIP' in row):
            continue
        sci_targets.append(row)

    sci_targets_table = ascii.read(sci_targets)
    for row in sci_targets_table:
        name = row['col1']
        ra_hms = f"%02i:%02i:%05.2f" % (row['col2'], row['col3'], row['col4'])
        dec_dms = f"%02i:%02i:%05.2f" % (row['col5'], row['col6'], row['col7'])
        print(name, ra_hms, dec_dms)

        tell_star = query_gemini_website(ra_hms=ra_hms, dec_dms=dec_dms, name=name)

        if tell_star == -99:
            continue
        else:
            hip_name = f"HIP{tell_star['HIP']}"
            ra_format = tell_star['RA(J2000)'].replace(':',' ')
            dec_format = tell_star['DEC(J2000)'].replace(':', ' ')
            comments = f"source_name = {name} Vmag={tell_star['Vmag']} dist={tell_star['Dist(deg)']}"
            with open(inp_filename,'a') as f:

                f.write(f"\n{hip_name.ljust(15)} {ra_format} {dec_format}0 2000.0  # {comments}")