#! /home/python3.7/bin/python
import json
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import ascii
import numpy as np
from astropy.time import Time
from datetime import datetime
import argparse
from astropy.table import Column
import matplotlib
import pandas as pd

matplotlib.use('Agg')


def init():
    matplotlib.rcParams['xtick.minor.size'] = 6
    matplotlib.rcParams['xtick.major.size'] = 6
    matplotlib.rcParams['ytick.major.size'] = 6
    matplotlib.rcParams['ytick.minor.size'] = 6
    matplotlib.rcParams['lines.linewidth'] = 1.5
    matplotlib.rcParams['axes.linewidth'] = 1.5
    matplotlib.rcParams['font.size'] = 16
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['xtick.major.width'] = 2.
    matplotlib.rcParams['ytick.major.width'] = 2.
    matplotlib.rcParams['ytick.direction'] = 'in'
    matplotlib.rcParams['xtick.direction'] = 'in'


def plot_fluxes(d, c='red', l='ZTF-g', ymn=-200, ymx=200, magticks=True):
    jd0 = Time(datetime.utcnow()).jd
    ax1.plot(d['jd,'] - jd0, np.array(d['forcediffimflux_uJy'], dtype=float), '.',
             color=c)
    ax1.errorbar(d['jd,'] - jd0, np.array(d['forcediffimflux_uJy'], dtype=float),
                 yerr=np.array(d['forcediffimfluxunc_uJy'], dtype=float), fmt='.',
                 color=c, label=l)  # ,t['col26']

    ax1.set_xlabel('JD - %.2f [today]' % (jd0))
    ax1.set_ylabel('Flux ($\mu$ Jy)')
    ax1.axhline(0, linestyle=':', c='black')
    ly, my = ax1.set_ylim(ymn, ymx)

    if magticks:
        ax2 = ax1.twinx()
        ax2.set_ylim(ly, my)
        mags = np.array([19, 19.5, 20, 20.5, 21])
        ax2.set_yticks(3631 * 10 ** (-0.4 * mags) * 1e6)
        ax2.set_yticklabels(mags)
        ax2.set_ylabel('mag')
    return plt


def plot_mags(d, c='red', l='ZTF-r', ymn=18, ymx=22):
    snr = np.array(d['forcediffimflux,'], dtype=float) / np.array(
        d['forcediffimfluxunc,'], dtype=float)
    snt = 3
    snu = 5
    detmask = (snr > snt)
    limmask = np.invert(detmask)

    ddet_mags = -2.5 * np.log10(
        np.array(d[detmask]['forcediffimflux_uJy'], dtype=float) * 1e-6 / 3631)
    ddet_magunc = 1.0857 * np.array(d[detmask]['forcediffimfluxunc_uJy'],
                                    dtype=float) / np.array(
        d[detmask]['forcediffimflux_uJy'], dtype=float)
    plt.plot(d[detmask]['jd,'] - jd0, ddet_mags, '.', color=c)
    plt.errorbar(d[detmask]['jd,'] - jd0, ddet_mags, ddet_magunc, fmt='.', color=c,
                 label=l)  # ,t['col26']

    dlim_mags = -2.5 * np.log10(
        snu * np.array(d[limmask]['forcediffimfluxunc_uJy'], dtype=float) * 1e-6 / 3631)
    plt.plot(d[limmask]['jd,'] - jd0, dlim_mags, 'v', color=c)
    plt.ylim(ymx, ymn)
    return plt


def plot_mags_atlas(atlaslcfile, ymn='18', ymx='22', jd0=0):
    atlaslc = ascii.read(atlaslcfile)
    c = atlaslc[atlaslc['F'] == 'c']
    o = atlaslc[atlaslc['F'] == 'o']
    c_dets = (c['uJy'] / c['duJy'] > 3)
    c_lims = (c['uJy'] / c['duJy'] < 3)
    o_dets = (o['uJy'] / o['duJy'] > 3)
    o_lims = (o['uJy'] / o['duJy'] < 3)

    print('ATLAS-c', len(c[c_dets]))
    print('ATLAS-o', len(o[o_dets]))
    plt.plot(c[c_dets]['##MJD'] + 2400000.5 - jd0,
             -2.5 * np.log10(c[c_dets]['uJy'] * 1e-6 / 3631), '.', c='cyan',
             label='ATLAS-c')
    plt.errorbar(c[c_dets]['##MJD'] + 2400000.5 - jd0,
                 -2.5 * np.log10(c[c_dets]['uJy'] * 1e-6 / 3631),
                 yerr=1.086 * c[c_dets]['duJy'] / c[c_dets]['uJy'], fmt='.', c='cyan')
    plt.plot(o[o_dets]['##MJD'] + 2400000.5 - jd0,
             -2.5 * np.log10(o[o_dets]['uJy'] * 1e-6 / 3631), '.', c='pink',
             label='ATLAS-o')
    plt.errorbar(o[o_dets]['##MJD'] + 2400000.5 - jd0,
                 -2.5 * np.log10(o[o_dets]['uJy'] * 1e-6 / 3631),
                 yerr=1.086 * o[o_dets]['duJy'] / o[o_dets]['uJy'], fmt='.', c='pink')

    plt.plot(c[c_lims]['##MJD'] + 2400000.5 - jd0,
             -2.5 * np.log10(5 * c[c_lims]['uJy'] * 1e-6 / 3631), 'v', c='cyan',
             label='ATLAS-o')
    plt.plot(o[o_lims]['##MJD'] + 2400000.5 - jd0,
             -2.5 * np.log10(5 * o[o_lims]['uJy'] * 1e-6 / 3631), 'v', c='pink',
             label='ATLAS-o')
    plt.ylim(ymx, ymn)
    return plt


def wavg(group, avg_name, avg2_name, weight_name):
    """ http://stackoverflow.com/questions/10951341/pandas-dataframe-aggregate-function-using-multiple-columns
    In rare instance, we may not have weights, so just return the mean. Customize this if your business case
    should return otherwise.
    """
    f = group[avg_name]
    d = group[avg2_name]
    w = (group[weight_name]) ** (-2)
    try:
        return (f * w).sum() / w.sum(), (d * w).sum() / w.sum(), (w.sum()) ** (-0.5)
    except ZeroDivisionError:
        return d.mean(), -1


def npwavg(group, avg_name, weight_name):
    # print group[avg_name].values, group[weight_name].values
    try:
        return np.average(group[avg_name], weights=1. / group[weight_name])
    except:
        return -1


def binavg(df, binsize):
    bfluxa, bmjda, bsigfluxa = [], [], []
    n = len(df)  # number of observations
    imin = 0
    imax = imin + binsize
    while imax <= n:
        # print df.iloc[imin:imax]
        bflux, bmjd, bsigflux = wavg(df.iloc[imin:imax], 'forcediffimflux_uJy', 'jd,',
                                     'forcediffimfluxunc_uJy')
        bfluxa.append(bflux)
        bmjda.append(bmjd)
        bsigfluxa.append(bsigflux)
        # print bmjd, bflux, bsigflux
        imin = imin + binsize
        imax = imax + binsize
    return np.array(bmjda), np.array(bfluxa), np.array(bsigfluxa)


def wavg_df(group, avg_name='forcediffimflux_uJy', avg2_name='jd,',
            weight_name='forcediffimfluxunc_uJy'):
    """ http://stackoverflow.com/questions/10951341/pandas-dataframe-aggregate-function-using-multiple-columns
    In rare instance, we may not have weights, so just return the mean. Customize this if your business case
    should return otherwise.
    """
    f = group[avg_name]
    d = group[avg2_name]
    w = (group[weight_name]) ** (-2)
    try:
        return pd.Series(
            {'bflux_uJy': (f * w).sum() / w.sum(), 'bjd': (d * w).sum() / w.sum(),
             'bfluxunc_uJy': (w.sum()) ** (-0.5)})
    except ZeroDivisionError:
        return d.mean(), -1


def plot_binned(d, binsize, c='red', l='ZTF-r', ymn=-50, ymx=100, magticks=True):
    bdf = pd.DataFrame.from_records(d, columns=d.colnames)
    bdjd, bdf, bdfunc = binavg(bdf, binsize)
    ax1.errorbar(bdjd - jd0, bdf, yerr=bdfunc, fmt='.', c=c,
                 label='%s, %sd' % (l, binsize))

    ax1.set_xlabel('JD - %.2f [today]' % (jd0))
    ax1.set_ylabel('Flux ($\mu$ Jy)')
    ax1.axhline(0, linestyle=':', c='black')
    ly, my = ax1.set_ylim(ymn, ymx)

    if magticks:
        ax2 = ax1.twinx()
        ax2.set_ylim(ly, my)
        mags = np.array([19, 19.5, 20, 20.5, 21])
        ax2.set_yticks(3631 * 10 ** (-0.4 * mags) * 1e6)
        ax2.set_yticklabels(mags)
        ax2.set_ylabel('mag')
    return plt


def plot_binned_mags_df(d, binsize, c='red', l='ZTF-r', ymn=18, ymx=22,
                        jd0=Time(datetime.utcnow()).jd, write=False):
    bdf = pd.DataFrame.from_records(d, columns=d.colnames)
    dbins = np.arange(int(bdf['jd,'].min()), int(bdf['jd,'].max()) + binsize, binsize)
    print(int(bdf['jd,'].min()), int(bdf['jd,'].max()), dbins)
    bdf_binned = bdf.groupby(np.digitize(bdf['jd,'], dbins))
    bdf_bin_flux = bdf_binned.apply(wavg_df)

    snt = 3
    snu = 5

    bdf_bin_flux['maglim'] = -2.5 * np.log10(
        bdf_bin_flux['bfluxunc_uJy'] * 5 * 1e-6 / 3631)
    if write:
        bdf_bin_flux.to_csv(f'fps_lcs/binned_{l}.csv')

    bdjd, bdf, bdfunc = bdf_bin_flux['bjd'], bdf_bin_flux['bflux_uJy'], bdf_bin_flux[
        'bfluxunc_uJy']
    snr = bdf / bdfunc

    detmask = (snr > snt)
    limmask = np.invert(detmask)

    ddet_mags = -2.5 * np.log10(np.array(bdf[detmask], dtype=float) * 1e-6 / 3631)
    ddet_magunc = 1.0857 * np.array(bdfunc[detmask], dtype=float) / np.array(
        bdf[detmask], dtype=float)

    # print('Max at',bdjd[detmask].iloc[np.argmin(ddet_mags)],ddet_mags[np.argmin(ddet_mags)],np.argmin(ddet_mags))
    plt.plot(bdjd[detmask] - jd0, ddet_mags, '.', color=c)
    plt.errorbar(bdjd[detmask] - jd0, ddet_mags, ddet_magunc, fmt='.', color=c,
                 label=l)  # ,t['col26']

    dlim_mags = -2.5 * np.log10(
        snu * np.array(bdfunc[limmask], dtype=float) * 1e-6 / 3631)
    plt.plot(bdjd[limmask] - jd0, dlim_mags, 'v', color=c, alpha=0.4)
    plt.ylim(ymx, ymn)

    # print(bdjd[detmask]-jd0)
    # print(bdjd[limmask]-jd0)
    return plt


def plot_binned_mags(d, binsize, c='red', l='ZTF-r', ymn=18, ymx=22, magticks=True):
    bdf = pd.DataFrame.from_records(d, columns=d.colnames)
    bdjd, bdf, bdfunc = binavg(bdf, binsize)
    snr = bdf / bdfunc
    snt = 3
    snu = 5
    detmask = (snr > snt)
    limmask = np.invert(detmask)

    ddet_mags = -2.5 * np.log10(np.array(bdf[detmask], dtype=float) * 1e-6 / 3631)
    ddet_magunc = 1.0857 * np.array(bdfunc[detmask], dtype=float) / np.array(
        bdf[detmask], dtype=float)
    plt.plot(bdjd[detmask] - jd0, ddet_mags, '.', color=c)
    plt.errorbar(bdjd[detmask] - jd0, ddet_mags, ddet_magunc, fmt='.', color=c,
                 label=l)  # ,t['col26']

    dlim_mags = -2.5 * np.log10(
        snu * np.array(bdfunc[limmask], dtype=float) * 1e-6 / 3631)
    plt.plot(bdjd[limmask] - jd0, dlim_mags, 'v', color=c)
    plt.ylim(ymx, ymn)
    return plt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--ztf_file", type=str, help="ZTF lightcurve file")
    parser.add_argument("--pgir_file", type=str, help="PGIR lightcurve file")
    parser.add_argument("--atlas_file", type=str, help="ATLAS lightcurve file")
    parser.add_argument("--mag", action='store_true', default=False,
                        help="make only magnitudes as well")
    parser.add_argument("--bin", action='store_true', default=False,
                        help="Bin photometry")
    parser.add_argument("--rbin", type=int, default=1,
                        help="Binsize in days for r-band")
    parser.add_argument("--gbin", type=int, default=1,
                        help="Binsize in days for g-band")
    parser.add_argument("--ibin", type=int, default=1,
                        help="Binsize in days for i-band")
    parser.add_argument("--cbin", type=int, default=1,
                        help="Binsize in days for cyan-band")
    parser.add_argument("--obin", type=int, default=1,
                        help="Binsize in days for orange-band")
    parser.add_argument("--xmin", default=None, type=float, help="set xlim to axis")
    parser.add_argument("--xmax", default=1, type=float, help="set xlim to axis")
    parser.add_argument("--ymin", default=18, type=float, help="set ylim to axis")
    parser.add_argument("--ymax", default=22, type=float, help="set ylim to axis")
    parser.add_argument("--time_isot_vline", default=None, type=str,
                        help="ISOT time to plot vertical line at")
    parser.add_argument("--mask", action="store_true",
                        help="skip applying mask on images")

    args = parser.parse_args()
    jd0 = Time(datetime.utcnow()).jd
    mjd_vline = None
    if args.time_isot_vline is not None:
        mjd_vline = Time(args.time_isot_vline).mjd
    init()

    if args.ztf_file is None and args.pgir_file is None and args.atlas_file is None:
        print('Please enter at least one lightcurve file')
        exit(0)

    fig, ax1 = plt.subplots(figsize=(10, 7))
    if args.ztf_file is not None:
        ztflc = ascii.read(args.ztf_file, format='basic')
        print('Read %s' % (len(ztflc)))

        try:
            np.array(ztflc['forcediffimflux,'], dtype=float)
        except ValueError:
            ztflc = ztflc[ztflc['forcediffimflux,'] != 'null']

        if args.mask:
            ztflc = ztflc[(ztflc['infobitssci,'] == 0) & (ztflc['scisigpix,'] < 25) & (
                    ztflc['sciinpseeing,'] < 4)]
        print('Found %s good points.' % (len(ztflc)))

        ztflc.add_column(Column(name='forcediffimflux_uJy',
                                data=np.array(ztflc['forcediffimflux,'],
                                              dtype=float) * 10 ** (-0.4 * (np.array(
                                    ztflc['zpmaginpsci,']) + 48.6)) / 1e-29))
        ztflc.add_column(Column(name='forcediffimfluxunc_uJy',
                                data=np.array(ztflc['forcediffimfluxunc,'],
                                              dtype=float) * 10 ** (-0.4 * (np.array(
                                    ztflc['zpmaginpsci,']) + 48.6)) / 1e-29))

        i = ztflc[ztflc['filter,'] == 'ZTF_i']
        r = ztflc[ztflc['filter,'] == 'ZTF_r']
        g = ztflc[ztflc['filter,'] == 'ZTF_g']
        jd0 = Time(datetime.utcnow()).jd

        if args.mag:
            if not args.bin:
                plt.figure(figsize=(10, 7))
                plt = plot_mags(r, c='red', l='ZTF-r', ymn=args.ymin, ymx=args.ymax)
                plt = plot_mags(i, c='brown', l='ZTF-i', ymn=args.ymin, ymx=args.ymax)
                plt = plot_mags(g, c='green', l='ZTF-g', ymn=args.ymin, ymx=args.ymax)
                # plt.axvline(specjd-jd0,linestyle='--',c='black')
                objname = args.ztf_file.split('.dat')[0]
                plt.legend(fontsize=12)
                plt.xlabel('JD - %i [today]' % (jd0))
                plt.ylabel('mag')
                if args.xmin is not None:
                    plt.xlim(args.xmin, args.xmax)

                if args.atlas_file is not None:
                    print('In Atlas')
                    atlaslcfile = args.atlas_file
                    plt = plot_mags_atlas(atlaslcfile, ymn=args.ymin, ymx=args.ymax,
                                          jd0=jd0)

                if mjd_vline is not None:
                    plt.axvline(mjd_vline + 2400000.5 - jd0,
                                linestyle='--', c='black', alpha=0.2)

                plt.savefig('%s_mags.pdf' % (objname), bbox_inches='tight')

            if args.bin:
                plt.figure(figsize=(10, 7))
                if len(r) > 0:
                    plt = plot_binned_mags_df(r, args.rbin, c='red', l='ZTF-r',
                                              ymn=args.ymin,
                                              ymx=args.ymax, jd0=jd0)
                if len(i) > 0:
                    plt = plot_binned_mags_df(i, args.ibin, c='brown', l='ZTF-i',
                                              ymn=args.ymin,
                                              ymx=args.ymax, jd0=jd0)
                if len(g) > 0:
                    plt = plot_binned_mags_df(g, args.gbin, c='green', l='ZTF-g',
                                              ymn=args.ymin,
                                              ymx=args.ymax, jd0=jd0)
                # plt.axvline(specjd-jd0,linestyle='--',c='black')
                objname = args.ztf_file.split('.dat')[0]
                plt.legend(fontsize=12)
                plt.xlabel('JD - %i [today]' % (jd0))
                plt.ylabel('mag')
                if args.xmin is not None:
                    plt.xlim(args.xmin, args.xmax)

                if args.atlas_file is not None:
                    print('In Atlas bin')
                    atlaslcfile = args.atlas_file
                    # objname = args.atlas_file.split('.dat')[0]

                    atlaslc = ascii.read(atlaslcfile)
                    atlaslc['jd,'] = atlaslc['##MJD'] + 2400000.5
                    atlaslc['forcediffimflux_uJy'] = np.array(atlaslc['uJy'],
                                                              dtype=float)
                    atlaslc['forcediffimfluxunc_uJy'] = np.array(atlaslc['duJy'],
                                                                 dtype=float)

                    c = atlaslc[atlaslc['F'] == 'c']
                    o = atlaslc[atlaslc['F'] == 'o']
                    plt = plot_binned_mags_df(c, args.cbin, c='cyan', l='ATLAS-c',
                                              ymn=args.ymin,
                                              ymx=args.ymax, jd0=jd0)
                    plt = plot_binned_mags_df(o, args.obin, c='orange', l='ATLAS-o',
                                              ymn=args.ymin,
                                              ymx=args.ymax, jd0=jd0)

                plt.legend()
                if mjd_vline is not None:
                    plt.axvline(mjd_vline + 2400000.5 - jd0,
                                linestyle='--', c='black', alpha=0.2)
                plt.savefig('%s_binned_mags.pdf' % (objname), bbox_inches='tight')

        plt = plot_fluxes(r, c='red', l='ZTF-r', magticks=False)
        plt = plot_fluxes(i, c='brown', l='ZTF-i', magticks=False)
        plt = plot_fluxes(g, c='green', l='ZTF-g', magticks=True)
        objname = args.ztf_file.split('.dat')[0]
        ax1.legend(fontsize=12)

        if args.bin:
            fig, ax1 = plt.subplots(figsize=(10, 7))
            plt = plot_binned(r, args.rbin, c='red', l='ZTF-r', magticks=False)
            plt = plot_binned(g, args.gbin, c='green', l='ZTF-g', magticks=False)
            plt = plot_binned(i, args.ibin, c='brown', l='ZTF-i', magticks=True)
            ax1.legend(fontsize=12)
            if mjd_vline is not None:
                plt.axvline(mjd_vline + 2400000.5 - jd0,
                            linestyle='--', c='black', alpha=0.2)
            plt.savefig('%s_binned.pdf' % (objname), bbox_inches='tight')

    if args.pgir_file is not None:
        pgirlcfile = args.pgir_file
        with open(pgirlcfile, 'r') as json_file:
            lc = json.load(json_file)

        jd0 = Time(datetime.utcnow()).jd
        for key in lc.keys():
            plt.plot(lc[key]['jd'] - jd0, lc[key]['psfflux'], '.', c='purple')
            plt.errorbar(lc[key]['jd'] - jd0, lc[key]['psfflux'],
                         lc[key]['psfflux_unc'], fmt='.', c='purple', label='PGIR J')

        objname = pgirlcfile.split('lc')[1].split('_')[1]

    if args.atlas_file is not None:
        atlaslcfile = args.atlas_file
        objname = args.atlas_file.split('.dat')[0]

        atlaslc = ascii.read(atlaslcfile)
        atlaslc['jd,'] = atlaslc['##MJD'] + 2400000.5
        atlaslc['forcediffimflux_uJy'] = np.array(atlaslc['uJy'], dtype=float)
        atlaslc['forcediffimfluxunc_uJy'] = np.array(atlaslc['duJy'], dtype=float)

        c = atlaslc[atlaslc['F'] == 'c']
        o = atlaslc[atlaslc['F'] == 'o']

        if not args.bin:
            ax1.plot(c['##MJD'] + 2400000.5 - jd0, c['uJy'], '.', c='cyan',
                     label='ATLAS-c')
            ax1.errorbar(c['##MJD'] + 2400000.5 - jd0, c['uJy'], c['duJy'], c='cyan',
                         linestyle='', fmt='.')

            ax1.plot(o['##MJD'] + 2400000.5 - jd0, o['uJy'], '.', c='pink',
                     label='ATLAS-o')
            ax1.errorbar(o['##MJD'] + 2400000.5 - jd0, o['uJy'], o['duJy'], c='pink',
                         linestyle='', fmt='.')

            ax1.legend(fontsize=9)

        if args.bin:
            fig, ax1 = plt.subplots(figsize=(10, 7))
            plt = plot_binned_mags_df(c, args.cbin, c='blue', l='ATLAS-c',
                                      ymn=args.ymin, ymx=args.ymax,
                                      jd0=jd0)
            plt = plot_binned_mags_df(o, args.obin, c='orange', l='ATLAS-o',
                                      ymn=args.ymin, ymx=args.ymax,
                                      jd0=jd0)
            if mjd_vline is not None:
                plt.axvline(mjd_vline + 2400000.5 - jd0,
                            linestyle='--', c='black', alpha=0.2)
            if args.xmin is not None:
                plt.xlim(args.xmin, args.xmax)
            plt.savefig('%s_binned_mags.pdf' % (objname), bbox_inches='tight')

    if args.xmin is not None:
        plt.xlim(args.xmin, args.xmax)

    if mjd_vline is not None:
        plt.axvline(mjd_vline + 2400000.5 - jd0,
                    linestyle='--', c='black', alpha=0.2)
    plt.savefig(f'{objname}.pdf', bbox_inches='tight')
