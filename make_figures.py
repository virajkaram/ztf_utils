import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from glob import glob
import numpy as np
from matplotlib import gridspec
import json
from astropy.time import Time


def init():
    matplotlib.rcParams['font.family']= 'sans-serif'
    #matplotlib.rcParams['xtick.major.width']= 2.
    #matplotlib.rcParams['ytick.major.width']= 2.
    matplotlib.rcParams['ytick.direction']='in'
    matplotlib.rcParams['xtick.direction']='in'

init()

def make_figure_lcs(lcfile,specfile,speclines=None,mol_bands=None,jd0=None,savedir='lc_pngs',uplim=None,lolim=None):
    fig = plt.figure(figsize=(10,6))

    with open(lcfile,'r') as json_file:
        lc = json.load(json_file)
    
    spec = fits.open(specfile)
    specdate = Time(spec[0].header['UTSHUT']).jd
    gs = gridspec.GridSpec(1,1)
    ax1 = fig.add_subplot(gs[0])

    if jd0 is None:
        jd0 = lc['epnew0']['jd']

    for key in lc.keys():
        if lc[key]['mag_ap3'] != -99 and lc[key]['magsnr_ap3']>3:
            ax1.errorbar(lc[key]['jd']-jd0,lc[key]['mag_ap3'],yerr=lc[key]['magerr_ap3'],c='purple',fmt='.',markersize=20)
        
        if lc[key]['maglim_ap3'] !=-99 and lc[key]['magsnr_ap3']<3:
            ax1.plot(lc[key]['jd']-jd0,lc[key]['maglim_ap3'],'v',c='purple')
            
    plt.axvline(specdate-jd0,color='black')
    plt.xlabel(r'JD $-$ %s'%(jd0),size=15)        
    plt.ylabel(r'mag',size=15)
    plt.gca().invert_yaxis()
    ax1.tick_params(size=15,labelsize=15)
    plt.tight_layout()
    plt.savefig(r'%s/lc_%s.png'%(savedir,specfile.split('/')[-1].split('.fits')[0]))
    plt.close()



def make_figure(lcfile,specfile,speclines=None,mol_bands=None,jd0=None,savedir='figures/He_zoom',uplim=None,lolim=None):
    fig = plt.figure(figsize=(10,6))

    with open(lcfile,'r') as json_file:
        lc = json.load(json_file)
    
    spec = fits.open(specfile)
    specdate = Time(spec[0].header['UTSHUT']).jd
    gs = gridspec.GridSpec(2,1)
    ax1 = fig.add_subplot(gs[0])

    if jd0 is None:
        jd0 = lc['epnew0']['jd']

    for key in lc.keys():
        if lc[key]['mag_ap3'] != -99 and lc[key]['magsnr_ap3']>3:
            ax1.errorbar(lc[key]['jd']-jd0,lc[key]['mag_ap3'],yerr=lc[key]['magerr_ap3'],c='purple',fmt='.',markersize=20)
        
        if lc[key]['maglim_ap3'] !=-99 and lc[key]['magsnr_ap3']<3:
            ax1.plot(lc[key]['jd']-jd0,lc[key]['maglim_ap3'],'v',c='purple')
            
    plt.axvline(specdate-jd0,color='black')
    plt.xlabel(r'JD $-$ %s'%(jd0),size=15)        
    plt.ylabel(r'mag',size=15)
    plt.gca().invert_yaxis()
    ax1.tick_params(size=15,labelsize=15)

    spec = fits.open(specfile)
    data = spec[0].data
    wavs = np.ndarray.flatten(np.array([data[3][0],data[2][0],data[1][0],data[0][0]]))
    fluxes = np.ndarray.flatten(np.array([data[3][1],data[2][1],data[1][1],data[0][1]]))    
    wavmask = ((wavs<1.46) & (wavs>1.35)) | ((wavs<1.93) & (wavs>1.8)) 
    wavmask = np.invert(wavmask)
    #cons = cons #+ 5e-14
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(wavs[wavmask],fluxes[wavmask],linewidth=0.7,c='r')


    plt.minorticks_on()
    if uplim is not None and lolim is not None:
        plt.ylim(lolim,uplim)
    plt.xlabel(r'$\lambda$',size=10)
    plt.ylabel(r'$f_{\lambda}$',size=10)

    ymn,ymx = plt.ylim()
    if speclines is not None:
        #ax2.vlines(speclines[np.where(speclines['Element']=='Fei')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='black',linewidth=0.5,alpha=0.3)
        ax2.vlines(speclines[np.where(speclines['Element']=='Ci')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='g',linewidth=0.5,alpha=0.3)
        ax2.vlines(H_lines,ymin=ymn,ymax=ymx,linestyle='--',color='brown',linewidth=0.5,alpha=0.7,label='HI')
        ax2.vlines(HeI,ymin=ymn,ymax=ymx,linestyle='--',color='black',linewidth=0.7,alpha=0.8)
        #ax2.vlines(speclines[np.where(speclines['Element']=='Cai')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='g',linewidth=0.5,alpha=0.3)
        #ax2.vlines(speclines[np.where(speclines['Element']=='Nai')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='orange',linewidth=0.5,alpha=0.7,label='Nai')
        #ax2.vlines(speclines[np.where(speclines['Element']=='Sii')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='yellow',linewidth=0.5,alpha=0.7,label='Nai')        

    if mol_bands is not None:
        #ax2.vlines(C12O16,ymin=ymn,ymax=ymx,linestyle='--',color='blue',linewidth=0.5,alpha=0.7,label='C12O16')
        #ax2.vlines(C12O18,ymin=ymn,ymax=ymx,linestyle='--',color='cyan',linewidth=0.5,alpha=0.7,label='C12O18')
        ax2.vlines(mol_bands[np.where(mol_bands['Molecule']=='CN')]['wav'],ymin=ymn,ymax=ymx,linestyle='--',color='b',linewidth=0.5,alpha=0.3)
        #ax2.vlines(mol_bands[np.where(mol_bands['Molecule']=='C2')]['wav'],ymin=ymn,ymax=ymx,linestyle='--',color='b',linewidth=0.5,alpha=0.7,label='C2')

    ax2.tick_params(size=15,labelsize=15)
    plt.xlim(1.03,1.12)
    plt.savefig(r'%s/lc_%s.pdf'%(savedir,specfile.split('/')[-1].split('.fits')[0]))
    plt.close()



def make_zoom_figure(specfile,speclines=None,mol_bands=None,jd0=None,savedir='proc_irtf/zoom',uplim=None,lolim=None):

    fig = plt.figure(figsize=(8,6))
    '''
    with open(lcfile,'r') as json_file:
        lc = json.load(json_file)
    
    spec = fits.open(specfile)
    specdate = Time(spec[0].header['UTSHUT']).jd
    '''
    gs = gridspec.GridSpec(4,1)
    ax1 = fig.add_subplot(gs[0])
    '''
    if jd0 is None:
        jd0 = lc['epnew0']['jd']

    for key in lc.keys():
        if lc[key]['mag_ap3'] != -99 and lc[key]['magsnr_ap3']>3:
            ax1.errorbar(lc[key]['jd']-jd0,lc[key]['mag_ap3'],yerr=lc[key]['magerr_ap3'],c='purple',fmt='.',markersize=10)
        
        if lc[key]['maglim_ap3'] !=-99 and lc[key]['magsnr_ap3']<3:
            ax1.plot(lc[key]['jd']-jd0,lc[key]['maglim_ap3'],'v',c='purple')
            
    plt.axvline(specdate-jd0,color='black')
    plt.xlabel(r'JD $-$ %s'%(jd0),size=5)        
    plt.ylabel(r'mag',size=5)
    plt.gca().invert_yaxis()
    ax1.tick_params(size=2,labelsize=8)
    '''
    spec = fits.open(specfile)
    data = spec[0].data
    wavs = np.ndarray.flatten(np.array([data[3][0],data[2][0],data[1][0],data[0][0]]))
    fluxes = np.ndarray.flatten(np.array([data[3][1],data[2][1],data[1][1],data[0][1]]))    
    wavmask = ((wavs<1.12) & (wavs>1.03)) #| ((wavs<1.93) & (wavs>1.8)) 
    #wavmask = np.invert(wavmask)
    #cons = cons #+ 5e-14
    data = np.vstack([[wavs*1e4, fluxes]]).T
    with open('%s.txt'%(sname),'w') as f:
        np.savetxt(f,data)
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(wavs[wavmask],fluxes[wavmask],linewidth=0.7,c='r')


    plt.minorticks_on()
    if uplim is not None and lolim is not None:
        plt.ylim(lolim,uplim)
    plt.xlabel(r'$\lambda$',size=5)
    plt.ylabel(r'$f_{\lambda}$',size=5)

    ymn,ymx = plt.ylim()
    if speclines is not None:
        #ax2.vlines(speclines[np.where(speclines['Element']=='Fei')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='black',linewidth=0.5,alpha=0.3)
        ax2.vlines(speclines[np.where(speclines['Element']=='Ci')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='g',linewidth=0.5,alpha=0.3)
        ax2.vlines(H_lines,ymin=ymn,ymax=ymx,linestyle='--',color='brown',linewidth=0.5,alpha=0.7,label='HI')
        ax2.vlines(HeI,ymin=ymn,ymax=ymx,linestyle='--',color='black',linewidth=0.7,alpha=0.8)
        #ax2.vlines(speclines[np.where(speclines['Element']=='Cai')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='g',linewidth=0.5,alpha=0.3)
        #ax2.vlines(speclines[np.where(speclines['Element']=='Nai')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='orange',linewidth=0.5,alpha=0.7,label='Nai')
        #ax2.vlines(speclines[np.where(speclines['Element']=='Sii')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='yellow',linewidth=0.5,alpha=0.7,label='Nai')        

    if mol_bands is not None:
        ax2.vlines(C12O16,ymin=ymn,ymax=ymx,linestyle='--',color='blue',linewidth=0.5,alpha=0.7,label='C12O16')
        ax2.vlines(C12O18,ymin=ymn,ymax=ymx,linestyle='--',color='cyan',linewidth=0.5,alpha=0.7,label='C12O18')
        ax2.vlines(mol_bands[np.where(mol_bands['Molecule']=='CN')]['wav'],ymin=ymn,ymax=ymx,linestyle='--',color='b',linewidth=0.5,alpha=0.3)
        #ax2.vlines(mol_bands[np.where(mol_bands['Molecule']=='C2')]['wav'],ymin=ymn,ymax=ymx,linestyle='--',color='b',linewidth=0.5,alpha=0.7,label='C2')

    ax2.tick_params(size=5,labelsize=8)
    plt.xlim(1.03,1.12)

    ax3 = fig.add_subplot(gs[2])
    wavmask = ((wavs<1.8) & (wavs>1.5))
    ax3.plot(wavs[wavmask],fluxes[wavmask],linewidth=0.7,c='r')


    plt.minorticks_on()
    if uplim is not None and lolim is not None:
        plt.ylim(lolim,uplim)
    plt.xlabel(r'$\lambda$',size=5)
    plt.ylabel(r'$f_{\lambda}$',size=5)

    ymn,ymx = plt.ylim()
    if speclines is not None:
        #ax3.vlines(speclines[np.where(speclines['Element']=='Fei')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='black',linewidth=0.5,alpha=0.3)
        ax3.vlines(speclines[np.where(speclines['Element']=='Ci')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='g',linewidth=0.5,alpha=0.3)
        ax3.vlines(H_lines,ymin=ymn,ymax=ymx,linestyle='--',color='brown',linewidth=0.5,alpha=0.7,label='HI')
        ax3.vlines(HeI,ymin=ymn,ymax=ymx,linestyle='--',color='black',linewidth=0.7,alpha=0.8)
        #ax3.vlines(speclines[np.where(speclines['Element']=='Cai')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='g',linewidth=0.5,alpha=0.3)
        #ax3.vlines(speclines[np.where(speclines['Element']=='Nai')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='orange',linewidth=0.5,alpha=0.7,label='Nai')
        #ax3.vlines(speclines[np.where(speclines['Element']=='Sii')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='yellow',linewidth=0.5,alpha=0.7,label='Nai')        

    if mol_bands is not None:
        ax3.vlines(C12O16,ymin=ymn,ymax=ymx,linestyle='--',color='blue',linewidth=0.5,alpha=0.7,label='C12O16')
        ax3.vlines(C12O18,ymin=ymn,ymax=ymx,linestyle='--',color='cyan',linewidth=0.5,alpha=0.7,label='C12O18')
        ax3.vlines(mol_bands[np.where(mol_bands['Molecule']=='CN')]['wav'],ymin=ymn,ymax=ymx,linestyle='--',color='b',linewidth=0.5,alpha=0.3)
        #ax3.vlines(mol_bands[np.where(mol_bands['Molecule']=='C2')]['wav'],ymin=ymn,ymax=ymx,linestyle='--',color='b',linewidth=0.5,alpha=0.7,label='C2')

    ax3.tick_params(size=5,labelsize=8)
    plt.xlim(1.5,1.8)


    ax4 = fig.add_subplot(gs[3])
    wavmask = ((wavs<2.43) & (wavs>2.12))
    ax4.plot(wavs[wavmask],fluxes[wavmask],linewidth=0.7,c='r')


    plt.minorticks_on()
    if uplim is not None and lolim is not None:
        plt.ylim(lolim,uplim)
    plt.xlabel(r'$\lambda$',size=5)
    plt.ylabel(r'$f_{\lambda}$',size=5)

    ymn,ymx = plt.ylim()
    if speclines is not None:
        #ax4.vlines(speclines[np.where(speclines['Element']=='Fei')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='black',linewidth=0.5,alpha=0.3)
        ax4.vlines(speclines[np.where(speclines['Element']=='Ci')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='g',linewidth=0.5,alpha=0.3)
        ax4.vlines(H_lines,ymin=ymn,ymax=ymx,linestyle='--',color='brown',linewidth=0.5,alpha=0.7,label='HI')
        ax4.vlines(HeI,ymin=ymn,ymax=ymx,linestyle='--',color='black',linewidth=0.7,alpha=0.8)
        #ax4.vlines(speclines[np.where(speclines['Element']=='Cai')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='g',linewidth=0.5,alpha=0.3)
        #ax4.vlines(speclines[np.where(speclines['Element']=='Nai')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='orange',linewidth=0.5,alpha=0.7,label='Nai')
        #ax4.vlines(speclines[np.where(speclines['Element']=='Sii')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='yellow',linewidth=0.5,alpha=0.7,label='Nai')        

    if mol_bands is not None:
        ax4.vlines(C12O16,ymin=ymn,ymax=ymx,linestyle='--',color='blue',linewidth=0.5,alpha=0.7,label='C12O16')
        ax4.vlines(C12O18,ymin=ymn,ymax=ymx,linestyle='--',color='cyan',linewidth=0.5,alpha=0.7,label='C12O18')
        ax4.vlines(mol_bands[np.where(mol_bands['Molecule']=='CN')]['wav'],ymin=ymn,ymax=ymx,linestyle='--',color='b',linewidth=0.5,alpha=0.3)
        #ax4.vlines(mol_bands[np.where(mol_bands['Molecule']=='C2')]['wav'],ymin=ymn,ymax=ymx,linestyle='--',color='b',linewidth=0.5,alpha=0.7,label='C2')

    ax4.tick_params(size=5,labelsize=8)
    plt.xlim(2.12,2.43)

    plt.savefig(r'%s/lc_%s.pdf'%(savedir,specfile.split('/')[-1].split('.fits')[0]))
    plt.close()


def make_zoom_figure_merged(lcfile,specfile,speclines=None,mol_bands=None,jd0=None,savedir='proc_irtf/zoom',uplim=None,lolim=None):
    fig = plt.figure(figsize=(8,6))

    with open(lcfile,'r') as json_file:
        lc = json.load(json_file)
    
    spec = fits.open(specfile)
    specdate = Time(spec[0].header['DATE_OBS']).jd
    gs = gridspec.GridSpec(4,1)
    ax1 = fig.add_subplot(gs[0])

    if jd0 is None:
        jd0 = lc['epnew0']['jd']

    for key in lc.keys():
        if lc[key]['mag_ap3'] != -99 and lc[key]['magsnr_ap3']>3:
            ax1.errorbar(lc[key]['jd']-jd0,lc[key]['mag_ap3'],yerr=lc[key]['magerr_ap3'],c='purple',fmt='.',markersize=10)
        
        if lc[key]['maglim_ap3'] !=-99 and lc[key]['magsnr_ap3']<3:
            ax1.plot(lc[key]['jd']-jd0,lc[key]['maglim_ap3'],'v',c='purple')
            
    plt.axvline(specdate-jd0,color='black')
    plt.xlabel(r'JD $-$ %s'%(jd0),size=5)        
    plt.ylabel(r'mag',size=5)
    plt.gca().invert_yaxis()
    ax1.tick_params(size=2,labelsize=8)

    spec = fits.open(specfile)
    data = spec[0].data
    wavs = data[0]#np.ndarray.flatten(np.array([data[3][0],data[2][0],data[1][0],data[0][0]]))
    fluxes = data[1]#np.ndarray.flatten(np.array([data[3][1],data[2][1],data[1][1],data[0][1]]))    
    wavmask = ((wavs<1.12))# & (wavs>1.03)) #| ((wavs<1.93) & (wavs>1.8)) 
    #wavmask = np.invert(wavmask)
    #cons = cons #+ 5e-14
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(wavs[wavmask],fluxes[wavmask],linewidth=0.7,c='r')


    plt.minorticks_on()
    if uplim is not None and lolim is not None:
        plt.ylim(lolim,uplim)
    plt.xlabel(r'$\lambda$',size=5)
    plt.ylabel(r'$f_{\lambda}$',size=5)

    ymn,ymx = plt.ylim()
    if speclines is not None:
        #ax2.vlines(speclines[np.where(speclines['Element']=='Fei')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='black',linewidth=0.5,alpha=0.3)
        ax2.vlines(speclines[np.where(speclines['Element']=='Ci')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='g',linewidth=0.5,alpha=0.3)
        ax2.vlines(H_lines,ymin=ymn,ymax=ymx,linestyle='--',color='brown',linewidth=0.5,alpha=0.7,label='HI')
        ax2.vlines(HeI,ymin=ymn,ymax=ymx,linestyle='--',color='black',linewidth=0.7,alpha=0.8)
        #ax2.vlines(speclines[np.where(speclines['Element']=='Cai')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='g',linewidth=0.5,alpha=0.3)
        #ax2.vlines(speclines[np.where(speclines['Element']=='Nai')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='orange',linewidth=0.5,alpha=0.7,label='Nai')
        #ax2.vlines(speclines[np.where(speclines['Element']=='Sii')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='yellow',linewidth=0.5,alpha=0.7,label='Nai')        

    if mol_bands is not None:
        ax2.vlines(C12O16,ymin=ymn,ymax=ymx,linestyle='--',color='blue',linewidth=0.5,alpha=0.7,label='C12O16')
        ax2.vlines(C12O18,ymin=ymn,ymax=ymx,linestyle='--',color='cyan',linewidth=0.5,alpha=0.7,label='C12O18')
        ax2.vlines(mol_bands[np.where(mol_bands['Molecule']=='CN')]['wav'],ymin=ymn,ymax=ymx,linestyle='--',color='b',linewidth=0.5,alpha=0.3)
        #ax2.vlines(mol_bands[np.where(mol_bands['Molecule']=='C2')]['wav'],ymin=ymn,ymax=ymx,linestyle='--',color='b',linewidth=0.5,alpha=0.7,label='C2')

    ax2.tick_params(size=5,labelsize=8)
    plt.xlim(wavs[0],1.12)

    ax3 = fig.add_subplot(gs[2])
    wavmask = ((wavs<1.8) & (wavs>1.5))
    ax3.plot(wavs[wavmask],fluxes[wavmask],linewidth=0.7,c='r')


    plt.minorticks_on()
    if uplim is not None and lolim is not None:
        plt.ylim(lolim,uplim)
    plt.xlabel(r'$\lambda$',size=5)
    plt.ylabel(r'$f_{\lambda}$',size=5)

    ymn,ymx = plt.ylim()
    if speclines is not None:
        #ax3.vlines(speclines[np.where(speclines['Element']=='Fei')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='black',linewidth=0.5,alpha=0.3)
        ax3.vlines(speclines[np.where(speclines['Element']=='Ci')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='g',linewidth=0.5,alpha=0.3)
        ax3.vlines(H_lines,ymin=ymn,ymax=ymx,linestyle='--',color='brown',linewidth=0.5,alpha=0.7,label='HI')
        ax3.vlines(HeI,ymin=ymn,ymax=ymx,linestyle='--',color='black',linewidth=0.7,alpha=0.8)
        #ax3.vlines(speclines[np.where(speclines['Element']=='Cai')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='g',linewidth=0.5,alpha=0.3)
        #ax3.vlines(speclines[np.where(speclines['Element']=='Nai')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='orange',linewidth=0.5,alpha=0.7,label='Nai')
        #ax3.vlines(speclines[np.where(speclines['Element']=='Sii')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='yellow',linewidth=0.5,alpha=0.7,label='Nai')        

    if mol_bands is not None:
        ax3.vlines(C12O16,ymin=ymn,ymax=ymx,linestyle='--',color='blue',linewidth=0.5,alpha=0.7,label='C12O16')
        ax3.vlines(C12O18,ymin=ymn,ymax=ymx,linestyle='--',color='cyan',linewidth=0.5,alpha=0.7,label='C12O18')
        ax3.vlines(mol_bands[np.where(mol_bands['Molecule']=='CN')]['wav'],ymin=ymn,ymax=ymx,linestyle='--',color='b',linewidth=0.5,alpha=0.3)
        #ax3.vlines(mol_bands[np.where(mol_bands['Molecule']=='C2')]['wav'],ymin=ymn,ymax=ymx,linestyle='--',color='b',linewidth=0.5,alpha=0.7,label='C2')

    ax3.tick_params(size=5,labelsize=8)
    plt.xlim(1.5,1.8)


    ax4 = fig.add_subplot(gs[3])
    wavmask = ((wavs<2.43) & (wavs>2.12))
    ax4.plot(wavs[wavmask],fluxes[wavmask],linewidth=0.7,c='r')


    plt.minorticks_on()
    if uplim is not None and lolim is not None:
        plt.ylim(lolim,uplim)
    plt.xlabel(r'$\lambda$',size=5)
    plt.ylabel(r'$f_{\lambda}$',size=5)

    ymn,ymx = plt.ylim()
    if speclines is not None:
        #ax4.vlines(speclines[np.where(speclines['Element']=='Fei')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='black',linewidth=0.5,alpha=0.3)
        ax4.vlines(speclines[np.where(speclines['Element']=='Ci')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='g',linewidth=0.5,alpha=0.3)
        ax4.vlines(H_lines,ymin=ymn,ymax=ymx,linestyle='--',color='brown',linewidth=0.5,alpha=0.7,label='HI')
        ax4.vlines(HeI,ymin=ymn,ymax=ymx,linestyle='--',color='black',linewidth=0.7,alpha=0.8)
        #ax4.vlines(speclines[np.where(speclines['Element']=='Cai')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='g',linewidth=0.5,alpha=0.3)
        #ax4.vlines(speclines[np.where(speclines['Element']=='Nai')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='orange',linewidth=0.5,alpha=0.7,label='Nai')
        #ax4.vlines(speclines[np.where(speclines['Element']=='Sii')]['Wav'],ymin=ymn,ymax=ymx,linestyle='--',color='yellow',linewidth=0.5,alpha=0.7,label='Nai')        

    if mol_bands is not None:
        ax4.vlines(C12O16,ymin=ymn,ymax=ymx,linestyle='--',color='blue',linewidth=0.5,alpha=0.7,label='C12O16')
        ax4.vlines(C12O18,ymin=ymn,ymax=ymx,linestyle='--',color='cyan',linewidth=0.5,alpha=0.7,label='C12O18')
        ax4.vlines(mol_bands[np.where(mol_bands['Molecule']=='CN')]['wav'],ymin=ymn,ymax=ymx,linestyle='--',color='b',linewidth=0.5,alpha=0.3)
        #ax4.vlines(mol_bands[np.where(mol_bands['Molecule']=='C2')]['wav'],ymin=ymn,ymax=ymx,linestyle='--',color='b',linewidth=0.5,alpha=0.7,label='C2')

    ax4.tick_params(size=5,labelsize=8)
    plt.xlim(2.12,2.43)

    plt.savefig(r'%s/lc_%s.pdf'%(savedir,specfile.split('/')[-1].split('.fits')[0]))
    plt.close()


H_lines = np.sort(np.array([0.410,0.434,0.486,0.656,1.875,1.282,1.094,1.005,0.954,2.625,2.166,1.944,1.817,1.458,2.279,1.94509,1.81791,1.73669,1.68111,1.64117,1.61137,1.58849,1.57050,1.55607,1.54431,1.53460,1.52647,1.51960,1.51374]))
C12O16  = np.array([2.2935,2.3227,2.3525,2.3829,2.4141])
C12O18 = np.array([2.3492,2.3783,2.4385])
C12O17 = np.array([2.3226,2.3517,2.3815,2.4119])
C13O16 = np.array([2.3448,2.3739,2.4037,2.4341])
HeI = np.array([1.083])#,1.701,1.278,2.059,2.113,2.165])


#speclines = ascii.read('sun_spectral_lines.csv')
#mol_bands = ascii.read('molecular_bands.csv')
'''
speclist = glob('proc/spec_rcb*')
#ids = [327,2822,2844,2874,23,267,2046,4018,4020]
for s in speclist:
    print(s)
    rcbid = s.split('_')[1].split('rcb')[-1]
    if s.split('_')[-1] == 'tellspec.fits':
        continue
    
    lcfile = glob('reqid_373/lc_%s_*.json'%(rcbid))[0]
    make_zoom_figure(lcfile,s,speclines=speclines,mol_bands=mol_bands)
    print("Done with %s light_curve %s"%(s,lcfile))
'''
'''
speclist = glob('proc/spec_known_rcb*')
for s in speclist:
    #print(s)
    rcbid = s.split('_')[-2].split('rcb')[-1]
    if s.split('_')[-1] == 'tellspec.fits':
        continue
    lcfile = glob('reqid_374/lc_%s_*.json'%(rcbid))[0]
    make_zoom_figure(lcfile,s,speclines=speclines,mol_bands=mol_bands)
    print("Done with %s light_curve %s"%(s,lcfile))
'''
'''
speclist = glob('proc/spec_rcb*')
for s in speclist:
    print(s)
    rcbid = s.split('_')[-2].split('rcb')[-1]
    if s.split('_')[-1] == 'tellspec.fits':
        continue
    lcfile = glob('reqid_358/lc_%s_*.json'%(rcbid))[0]
    make_zoom_figure(lcfile,s,speclines=speclines,mol_bands=mol_bands)
    print("Done with %s light_curve %s"%(s,lcfile))
'''
'''
speclist = glob('proc_irtf/spec_rcb*merged.fits')
print(speclist)
for s in speclist:
    rcbid = s.split('_')[-3].split('rcb')[-1]
    if s.split('_')[-1] == 'tellspec.fits':
        continue
    try:
        lcfile = glob('reqid_431/lc_%s_*.json'%(rcbid))[0]
    except IndexError:
        try:
            lcfile = glob('reqid_373/lc_%s_*.json'%(rcbid))[0]
        except:
            print('lc not found, skipping')
            continue
    make_zoom_figure_merged(lcfile,s,speclines=speclines,mol_bands=mol_bands)
    print("Done with %s light_curve %s"%(s,lcfile))
'''
sname = 'ztf20acdqjeq.fits'
make_zoom_figure(sname,savedir='/home/viraj/ZTF/')