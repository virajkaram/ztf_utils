#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import requests
from dustmaps.config import config
config['data_dir'] = '/home/viraj/dustmaps'
from dustmaps.sfd import SFDQuery
from astropy.coordinates import SkyCoord
import astropy.units as u
import argparse
from astropy.time import Time
from datetime import datetime
import json

# In[2]:


def gen_collc_interp(photometry,flux=False):
    if flux:
        data = np.array([(x['mjd'],x['flux'],x['filter'][-1],-2.5*np.log10(5*x['fluxerr']*1e-6/3631),x['fluxerr']) for x in photometry])
        data[:,4][data[:,1]!=None] = 2.5*np.array(data[:,4][data[:,1]!=None],dtype=float)/np.array(data[:,1][data[:,1]!=None],dtype=float)
        data[:,4][data[:,1]==None] = None
        data[:,1][data[:,1]!=None] = -2.5*np.log10(np.array(data[:,1][data[:,1]!=None],dtype=float)*1e-6/3631)

    else:    
        data = np.array([(x['mjd'],x['mag'],x['filter'][-1],x['limiting_mag'],x['magerr']) for x in photometry])
    data[:,4][data[:,4]==None] = 0

    gpoints = data[(data[:,2]=='g')]
    rpoints = data[(data[:,2]=='r')]
    gds = np.array(gpoints[:,0],dtype=float)
    rds = np.array(rpoints[:,0],dtype=float)
    rdetds = np.array(rds[(rpoints[:,1] != None)],dtype=float)
    gdetds = np.array(gds[(gpoints[:,1] != None)],dtype=float)
    
    if (len(gpoints)==0)|(len(rdetds)==0):
        return (np.array([0]),np.array([0]),np.array([0]),np.array([True]),np.array([True]),np.array([0]),np.array([0]))
    if len(rdetds)>len(gdetds):
        interpolate_points = rpoints
        fix_points = gpoints
    else:
        interpolate_points = gpoints
        fix_points = rpoints
    
    fix_ds = np.array(fix_points[:,0],dtype=float)
    fix_detds = np.array(fix_ds[(fix_points[:,1] != None)],dtype=float)
    fix_mags = np.array(fix_points[:,1],dtype=float)
    fix_magerrs = np.array(fix_points[:,4],dtype=float)
    fix_detmags = np.array(fix_points[:,1][(fix_points[:,1]!=None)],dtype=float)
    
    interpolate_ds = np.array(interpolate_points[:,0],dtype=float)
    interpolate_detds = np.array(interpolate_ds[(interpolate_points[:,1] != None)],dtype=float)
    interpolate_mags = np.array(interpolate_points[:,1],dtype=float)
    interpolate_magerrs = np.array(interpolate_points[:,4],dtype=float)
    interpolate_detmags = np.array(interpolate_points[:,1][(interpolate_points[:,1] != None)],dtype=float)
    
    interps=np.interp(fix_ds,interpolate_detds,interpolate_detmags)
    interps[fix_ds<np.min(interpolate_detds)] = np.nan

    clsinds = [np.argmin(np.abs(x-interpolate_ds)) for x in fix_ds]
    cldists = np.array([np.min(np.abs(x-interpolate_ds)) for x in fix_ds])
    rgs = interpolate_points[clsinds][:,1]
    rglims = np.array(interpolate_points[clsinds][:,3],dtype=float)
    ncnsismask = (interps<rglims)
    dtmask = (cldists<1)
    interp_limmask = (rgs==None)
    fix_limmask = (fix_points[:,1]==None)    
    interps[(ncnsismask&interp_limmask)] = rglims[(ncnsismask&interp_limmask)]
    interps[(dtmask&interp_limmask&fix_limmask)] = np.nan
        
    interperrs = interpolate_magerrs[clsinds]
        
    cols = fix_mags - interps
    cols[(fix_points[:,1]==None)] = np.array(fix_points[(fix_points[:,1]==None)][:,3],dtype=float) - interps[(fix_points[:,1]==None)] #- 2.5*np.log10(5./3) #seems like the lim_mags are already 5-sigma 
    colerrs = np.sqrt(interperrs**2 + fix_magerrs**2)
            
    fix_collims = (fix_points[:,1]== None) #is the g-band an upper limit? then g-r>g0-r
    interpolate_collims = interp_limmask&dtmask
    
    if len(rdetds)>len(gdetds):
        return (fix_ds,cols,colerrs,fix_collims,interpolate_collims,interpolate_detds,interpolate_detmags)
    
    else:
        return (fix_ds,-1*cols,colerrs,interpolate_collims,fix_collims,interpolate_detds,interpolate_detmags)


def is_subluminous(photometry,dm=None,flux=False):
    if flux:
        data = np.array([(x['mjd'],x['flux']) for x in photometry if x['flux'] !=None])
        data[:,1] = -2.5*np.log10(np.array(data[:,1],dtype=float)*1e-6/3631)
    else:    
        data = np.array([(x['mjd'],x['mag']) for x in photometry if x['mag']!=None])
    
    if dm==None:
        return True, data[:,0][np.argmin(data[:,1])]
    return np.min(data[:,1])-dm > -16, data[:,0][np.argmin(data[:,1])]
    

def lrn_filter(photometry,fore_gmr,dm=None,flux=False):
    redalways = False
    redlater = False
    subluminous = False
    bluepeak = False
    
    subluminous,maxday = is_subluminous(photometry,dm,flux=flux)
    d,c,cerr,gl,rl,dds,dmags = gen_collc_interp(photometry,flux=flux)
    c = c - fore_gmr
    #maxday = dds[np.argmin(dmags)]
    if len(c)==1 and np.isnan(c[0]):
        return False
    
    mdcdind = np.argmin(np.abs(d-maxday))
    maxdaycold = d[mdcdind]
    maxdaycol = c[mdcdind]
    '''
    redd = d[(c-cerr>0.7) & (d>maxdaycold) & (rl==False)]
    redlater = (len(redd)>1)
    if redlater:
        redlater = redlater & (np.min(redd)<maxdaycold+20)
    '''
    bluepeak = maxdaycol<0.5
    blued = d[(gl==False) & (c-cerr<0.5)]
    #print(maxday,blued)
    if len(blued) == 0:
        redlater = True
        redalways = True
    else:
        redlater = np.invert(np.any(blued>maxdaycold+20))
        
    #redalways = np.all(c[(rl==False)&(gl==False)]>0.7)    
    return (bluepeak&redlater | redalways) &subluminous


#groupids  41:rcf, 43:clu
def api(method, endpoint, data=None):
    with open('/Users/viraj/ztf_utils/secrets_fritz_token.json','r') as f:
        dat = json.load(f)
    token = dat['token']
    headers = {'Authorization': f'token {token}'}
    response = requests.request(method, endpoint, params=data, headers=headers)
    return response


if __name__ == '__main__':
    today = Time(datetime.utcnow()).mjd
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--start_date",type=float,default=58000,help="start date in MJD")
    parser.add_argument("--end_date",type=float,default=today,help="end date in MJD, default=today")
    parser.add_argument("--nslices",type=int,default=12,help="number of slices to query")
    parser.add_argument("--o",type=str,default='lrn_filter_%s.dat'%(today),help="output file to write results in")
    parser.add_argument("--prog_id",type=str,default='41',help="program id, example '41' or '41,43' for both")

    args = parser.parse_args()
    prog_id= args.prog_id
    nslices = args.nslices
    dates = np.linspace(args.start_date,args.end_date,nslices+1)
    dates = Time(dates,format='mjd')

    with open(args.o,'w') as f:
        f.write('Results for run from %s to %s'%(dates[0].isot,dates[-1].isot))

    for ind in range(len(dates)-1):
        start_date = dates[ind].isot
        end_date = dates[ind+1].isot
        data = {
            "group_ids" :prog_id,
            "startDate" : start_date,
            "endDate" : end_date,
            "includePhotometry":True
          }

        print('Querying sources in program %s, slice %i/%i from time %s to %s'%(prog_id,ind+1,args.nslices,start_date,end_date))
        response = api('GET', 'https://fritz.science/api/sources/',data=data)
        
        print(f'HTTP code: {response.status_code}, {response.reason}')
        print('Retrieved %s sources'%(len(response.json()['data']['sources'])))
        
        if response.status_code != 502:
            print('Timeout')
        if len(response.json()['data']['sources']) == 0:
            continue
        #Foreground extinction
        SF11 = {}
        SF11['g'] = 3.303
        SF11['r'] = 2.285
        SF11['i'] = 1.698
        
        crds = [(x['ra'],x['dec']) for x in response.json()['data']['sources']]
        crds = SkyCoord(crds,unit=(u.deg,u.deg))
        sfd = SFDQuery()
        ebvs = sfd(crds)
        
        
        Ag = (SF11['g'])*ebvs
        Ar = SF11['r']*ebvs
        Ai = SF11['i']*ebvs
        fore_gmr = Ag-Ar
        fore_gmi = Ag-Ai
        fore_rmi = Ar-Ai
        
        
        # In[112]:
        sources = response.json()['data']['sources']
        for k in range(len(sources)):
            source = sources[k]
            photometry = source['photometry']
            fgmr = fore_gmr[k]
            if len(photometry)==0:
                continue

            if lrn_filter(photometry,fgmr,dm =source['dm'],flux=True):
                print(source['id'])
                with open(args.o,'a') as f:
                    f.write('%s\n'%(source['id']))
        
        
        '''
        #Plot color curves
        for source in response.json()['data']['sources']:
            print(source['id'])
            photometry = source['photometry']
            d,c,cerr,gl,rl = gen_collc_interp(photometry,flux=True)
            plt.figure()
            plt.plot(d[gl],c[gl],'^')
            plt.plot(d[rl],c[rl],'v')
            plt.plot(d[np.invert(gl | rl)],c[np.invert(gl | rl)],'.')
            plt.ylabel('g - r')
            plt.xlabel('MJD')
            plt.title(source['id'])
            plt.savefig('col_lc/%s.pdf'%(source['id']))
            plt.close()
        
        '''
