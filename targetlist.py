from astropy.table import Table
import argparse
from astropy.coordinates import SkyCoord
import astropy.units as u


def make_palomar_targetlist(input_filename,name_colname='Toi_id',ra_colname='RA',dec_colname='Dec',comment_colname='Priority'):
    #palomar targetlist format

    file2 = Table.read(input_filename,format='csv')
    with open(input_filename.split('.')[0]+'_palomar.dat','w') as f:

        for i in range(len(file2)):
            coo = SkyCoord(ra=file2[i][ra_colname],dec=file2[i][dec_colname],unit=(u.deg,u.deg))
            rah, ram, ras = coo.ra.hms
            decd, decm, decs = coo.dec.dms
            try:
                name = 'RCB_%s'%(str(file2[i]["Toi_id"]))
            except:
                name = str(file2[i][name_colname])
                #name = 'RCB_%s'%(str(file2[i]['Name']))
            try:
                comment = '%s'%(file2[i]['J'])#['Tiss_web_id']/1000 +1)
            except:
                comment =''
            if coo.dec >=0:
                p = '+'
            else:
                p='-'
                decd = -1*decd
                decm = -1*decm
                decs = -1*decs
            
            try:
                priority = str(file2[i][comment_colname])
            except:
                priority = 'P4'
            f.write('%s%02i %02i %05.2f  %s%02i %02i %05.2f  2000.0  ! %s  J = %s\n'%(name.ljust(20),rah,ram,ras,p,decd,decm,decs,priority, comment))
       

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input_filename",type=str,help='Input targetlist in csv format')
    parser.add_argument("--name",type=str,default='Name',help='Column name for targetname')
    parser.add_argument("--ra",type=str,default='RA',help='Column name for RA')
    parser.add_argument("--dec",type=str,default='Dec',help='Column name for Dec')
    parser.add_argument("--comments",type=str,default=None,help='Column name for Dec')

    args = parser.parse_args()
    make_palomar_targetlist(args.input_filename, args.name, args.ra, args.dec, args.comments)
