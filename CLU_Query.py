# Author: Anastasios Tzanidakis
# Doc: Query all sources iN CLU using the GROWTH Marshall
import pickle
import numpy as np
from astropy.table import Table
from astropy.io import ascii
import os
import requests
import json
from datetime import datetime
from astropy.table import Table
from astropy.time import Time
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy import constants as const
from astropy.stats import *
from scipy import spatial
from scipy import interpolate
from scipy.interpolate import interp1d
from astroquery.irsa_dust import IrsaDust
import astropy.coordinates as coord
import json

import astropy.units as u
from scipy.optimize import curve_fit


# Username & Password Credentials
with open('/Users/viraj/ztf_utils/secrets_marshal.json','r') as f:
    dat = json.load(f)
user = dat['user']
passw = dat['pwd']

# What are the parameters we want to download from the GROWTH Marshal?
# ZTF_id, classification, number_of_spectra, {date_saved_iso}, z, ra, dec, peak_mag(red) or green when fail, d_CLU_galaxy, peak_abs_mag

ZTF_id, classification = [], [] # ZTF_id_# and classification
date_passed = [] # date passed the filter
n_spectra, z, ra, dec = [], [], [], [] # No. Spectra, distance and coords
peak_mag = [] # peak apparent magnitude in red
peak_abs_mag = [] # peak absolute magnitude in red
dist_to_galaxy = [] # distance from CLU galaxy
# Total: 9 variables (+1 extra peak abs_mag)

def linear_function(a, x, b):
    return (a*x+b)

# Program ID's you're part of: 0:CLU, 1:RCF
def what_programs_in(user, passw):
    """ This function will return the index and program you're in """
    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_programs.cgi', auth=(user, passw))
    programs = json.loads(r.text)
    programidx = -1
    print ("Programs you are a member of:")
    for index, program in enumerate(programs):
        print (program['programidx'], program['name'])

def dist_mod_mag(app_mag, distance):
    """ Calculate the distance modulus of cosmological sources """
    return (app_mag - 5*np.log10(distance)-25)

def redshift_to_distance(z):
    # Define cosmological constanst from Astropy
    cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
    return (cosmo.luminosity_distance(z))


prate = 0  # counting index for how many sources make it through the filter!

###### Begin extracting sources fro the GROWTH Marshall #######

programidx = 1 # Select program: RedTransients
r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_program_sources.cgi', auth=(user, passw), data = {'programidx': str(programidx)})
sources_clu = json.loads(r.text) # preliminary information on CLU objects
print(len(sources_clu))
massiveDict = {}
massiveDict['sources'] = []
for source in sources_clu[:2]:
    prate = prate+1
    s = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/source_summary.cgi', auth=(user, passw), data = {'sourceid': str(source['id'])})
    sourceDict = json.loads(s.text)
    massiveDict['sources'].append(sourceDict)
    print("Queried %s of %s"%(prate,len(sources_clu)))
    if prate%20 == 0:
        with open('redtransient_sources.json','w') as f:
            json.dump(massiveDict,f)

with open('redtransient_sources.json','w') as f:
    json.dump(massiveDict,f)

    '''
    spec = sourceDict['uploaded_spectra'] # information on spectroscopy
    spec_shape = np.shape(spec)
    annotation = sourceDict['annotations'] # acess the annotation dictionary
    auto_annotation = sourceDict['autoannotations'] # acess the annotation dictionary
    photometry_info = sourceDict['uploaded_photometry'] # information on photometry

    # Preliminary information --
    print ("   ")
    print ("   ")
    print ("   ")
    print ("###### Preliminary Information ######")
    print ("Currently examining %s"%source['name'])
    print ("Source URL:")
    print ("http://skipper.caltech.edu:8080/cgi-bin/growth/view_source.cgi?name=%s"%source['name'])
    print ("This source has N_spectra: %s"%spec_shape)
    print ("Source RA: %s"%sourceDict['ra'])
    print ("Source DEC %s"%sourceDict['dec'])
    ra_c, dec_c = sourceDict['ra'], sourceDict['dec'] # call ra, dec variables-- used for galactic exctinction corrections!
    print ("Source Saved Date: %s"%sourceDict["creationdate"]) # creation date!
    print ("__________________________")

    # Extract the distance from CLU_d_galaxy
    print ("###### Extracting CLU_d_galaxy_arcsec ######")
    temp_1 = [] # create temporary list
    for aa in auto_annotation:
        if aa['type']=='CLU_d_to_galaxy_arcsec':
            com = aa['comment'].split(" ") # split the output we find from the comments
            print ("Result of the galaxy d comment: %s"%com)
            print ("Distance to CLU Galaxy: %s"%com[0])
            temp_1.append(com[0])
    if np.shape(temp_1)[0]>0:
        print ("We actually have found a distance to a galaxy! wohooo! The CLU_d_galaxy we found: %s"%com[0])
        dist_to_galaxy.append(com[0])
    else:
        print ("Sorry couldn't find a distance to a galaxy!")
        dist_to_galaxy.append("None")
    print ("__________________________")


    # Extract redshift values
    print ("###### Extracting Redshift values ######")
    t_temp = [] # temporary list 1
    for ann in auto_annotation:
        if ann['type']=="CLU_z" or ann['type']=="CLU_Z" or ann['type']=="clu_z":
            comm_z = ann['comment'].split(" ")
            print ("We found something in the auto_annotation: %s"%comm_z[0])
            t_temp.append(comm_z[0])

    t_temp2 = [] # temporary list 2
    if np.shape(t_temp)[0]>=1: # if we found more than one in the comments!
        print ("Found value in the auto_annotation: %s"%t_temp[-1])
        z.append(t_temp[-1])

    elif np.shape(t_temp)[0]==0: # If we didn't find anything in the auto_annotation
        print ("No redshift was found... Now going to look in the comments")
        for ann2 in annotation:
            if ann2['type']=="redshift":
                z2 = ann2['comment']
                t_temp2.append(z2)
        if np.shape(t_temp2)[0]>=1:
            print (" We found something instead in the comments... this is the redshift I found: %s"%t_temp2[-1])
            z.append(t_temp2[-1])
        elif np.shape(t_temp2)[0]==0:
            if source['redshift']==None:
                print ("Sorry It seems I coudldn't find anything... Now appending a None")
                z.append("None")
            else:
                print ("Found something in the source header: %s"%source['redshift'])
                z.append(source['redshift'])
    print ("__________________________")

    # Extract saved dates from the GROWTH Marshall
    print ("###### Extracting Date-Saved values ######")
    t_creation = Time(sourceDict["creationdate"], format='iso')
    if t_creation==None:
        print ("Couldn't find a save date!")
        date_passed.append("None")
    else:
        print ("We have found a creation date! yay!")
        date_passed.append(t_creation.iso)
    print ("__________________________")

    # Extract Peak-appent magnitude sources from GROWTH Marshall

    print ("###### Extracting Peak-Apparent-Magnitudes values ######") # Fow now we will extrapolate peaks with percent growth!

    magnitudes = [photometry_info[i]['magpsf'] for i in range (0, len(photometry_info))] # append associated magnitudes
    mag_filter = [photometry_info[i]['filter'] for i in range (0, len(photometry_info))] # append assocaited magnitude filters
    #mag_errs = [photometry_info[i]['sigmamagpsf'] for i in range (0, len(photometry_info))] # erros will come into play at some point!
    jd_mag = [photometry_info[i]['jd'] for i in range (0, len(photometry_info))]
    isdiffpos = [photometry_info[i]['isdiffpos'] for i in range (0, len(photometry_info))]
    obsdate = [photometry_info[i]['obsdate'] for i in range (0, len(photometry_info))] # observed date!

    # create the magnitudes and mag_filters into numpy arrays
    magnitudes = np.asanyarray(magnitudes)
    mag_filter = np.asanyarray(mag_filter)
    #mag_errs = np.asanyarray(mag_errs) we will incorporate in the near future
    jd_mag = np.asanyarray(jd_mag)
    isdiffpos = np.asanyarray(isdiffpos)
    obsdate = np.asanyarray(obsdate)

    if len(magnitudes)==0 or len(obsdate)==0:
        # In the case where the magnitude dictionary is empty!
        print ("WARNING: There's nothing in this dictionary...")
        print ("Here's the N_source I found in those dict")
        print (len(magnitudes), len(obsdate), len(jd_mag))
        print ("Now going to append None to the peak_app_mag!")
        peak_mag.append("None")
        pass

    else:
        # If you find some data values in magnitudes let's move forward and select the r-band ones at isdifposs=True within reasonable values!

        # First let's select all red and green sources with mag <50(removing outliers) & are true detections (isdiffpos)
        red = np.where((mag_filter=='r') & (magnitudes>0) & (magnitudes<50) & (isdiffpos==True))

        # Okay now we have to decide if there's any red:
        N_red_sources = len(magnitudes[red])

        red_mag, red_jd = magnitudes[red], jd_mag[red] # assign to new variables!

        # Find the minima of the red_magnitude.
        minima_red = np.sort(red_mag)
        Number_of_minima_red = len(minima_red)

        if Number_of_minima_red > 0:
            try:
                minima_red = minima_red[1] # take the second hightest value!
            except:
                print ("Sorry we're gonna have to select the first index!")
                minima_red = minima_red[0]
        elif Number_of_minima_red <= 0: # if there's no minimum red...
            minima_red = False

        if minima_red==False:
            peak_mag.append("None")
            pass

        else:
            # Find the indext where the minimum red magnitude is!
            find_max = np.where(red_mag==minima_red)
            jdm = red_jd[find_max] # JD date where the minimum red is!
            jd_before_max = np.where(red_jd<jdm[0]) # find where is before peak
            jd_after_max = np.where(red_jd>jdm[0]) # find where after peak

            # Select before and after peak jd values!
            jd_bp, jd_ap = red_jd[jd_before_max]-jdm[0], red_jd[jd_after_max]-jdm[0] # note these are phases now!
            mr_bp, mr_ap = red_mag[jd_before_max], red_mag[jd_after_max]

            try:
                # Calculate a simple growth rate: max-min/min
                print (minima_red, np.max(mr_bp))
                growth = (minima_red-np.max(mr_bp))/np.max(mr_bp) # calculate growth rate
                gr = growth *100
                print ("Looks like your overal growth value was: %s"%gr)

                if gr < -2.5: # Generally -3% or less should be sufficient to catch all-if not all best cadence Lc's!
                    print ("Going into corrections now....")
                    skycoords = coord.SkyCoord(ra=ra_c*u.deg,dec=dec_c*u.deg, unit="deg") # covert to Astropy SkyCoord
                    table = IrsaDust.get_extinction_table(skycoords) # fetch dudst table
                    red_av = np.where(table['Filter_name']=='SDSS r') # select the SDSS Ar value
                    rav = (table['A_SFD'])[red_av] # Ar value
                    print ("Exiting corrections...")

                    print ("Gal magnitude correction!")
                    minima_red_gal_corr = minima_red -  rav.data # correct for Gal. Exctinciton using SDSS Av-r corr.
                    print ("Now appending peak magnitude: %s"%minima_red_gal_corr)
                    peak_mag.append(minima_red_gal_corr[0])
                else:
                    print ("Sorry your GR value is not sufficient enough to ensure that this is a true peak in the SNE. Now appending None")
                    peak_mag.append("None")
                    pass

            except:
                print ("Sorry, something went wrong with your LC extraction and couldn't calculate growth-rate. Now appending none...")
                peak_mag.append("None")
                pass


    # Extract Classifications from GROWTH Marshall
    print ("###### Extracting Source Classificaitons ######")

    print ("original header classification: %s"%sourceDict['classification'])
    ann_shape = np.shape(annotation)[0] # shape of annotation dictionary -- how manny annotations do we generally have

    # Run forloop to extrapolate classification and classification dates
    if ann_shape==1: # in the case there's only one annotation
        cl = sourceDict['classification'] # look at the header
        print ("No annotations posted... Looking for classification...")
        if cl==None:
            print ("No classification from the header... adding classification as None!")
            classification.append("None") # append None classification
        else:
            print ("We found something inside SourceDict['classification'] %s"%cl)
            classification.append(cl) # append that classification found in the header

    elif ann_shape>1:
        i = 0 # begin to keep an index
        types, types2, types3 = [], [], [] # types: classification temp list, types2: classification date temp list & jd and isot
        for ann in annotation: # loop through each annotation and find which one has the classification
            if ann['type']=="classification":
                print (ann['comment'])
                t = Time(ann['date'], format='isot') # record the time of each comment
                print (t)
                types.append(ann['comment']) # append to the temporary list
                if t==None: # check to make sure that the time frame is not empty!
                    print ("Empty time array!?")
                    types2.append("None")
                    types3.append("None")
                else: # if it's NOT empty then proceed to append those values!
                    types2.append(t.jd) # append jd time
                    types3.append(t.iso.split(" ")[0]) # append isot time
                i = i + 1 #add to the index to keep track of how many annotations had the classification
        if i>1: # if there's more than one classification
            print ("Found duplicate... We will select the most recent one : %s"%types[-1])
            classification.append(types[-1]) # return the msot recent one
        elif i==1: # if there's only one
            print ("Only one class was found... now appending :%s"%types)
            classification.append(types[-1])
        elif i==0: # if there are no classication
            header_class = sourceDict['classification']
            if header_class==None: # check to make sure the classification is not in the heade
                print ("No classification in comments nor header...")
                classification.append("None")
            else:
                print ("Alternatively we find something in the header of class: %s"%header_class)
                classification.append(header_class)
    else:
        print ("Didn't find anything! ")
        classification.append("None")
    print ("______________________")


    print ("Now appending values to emtpy lists...")
    # Append values to empty lists!
    ZTF_id.append(source['name']) # append name
    n_spectra.append(spec_shape[0]) # append number of spectra available
    ra.append(sourceDict['ra']) # append RA
    dec.append(sourceDict['dec']) # append DEC


    print ("##### Cecking all sources are the same length... ########")
    print (len(ZTF_id), len(n_spectra), len(ra), len(dec), len(z), len(classification), len(dist_to_galaxy), len(date_passed), len(peak_mag))


print ("Query is Complete!")

# Store everything in backup_data
ZTF_id, ra, dec, z, classification = np.asanyarray(ZTF_id), np.asanyarray(ra), np.asanyarray(dec), np.asanyarray(z), np.asanyarray(classification)
dist_to_galaxy, n_spectra, date_passed, peak_mag = np.asanyarray(dist_to_galaxy), np.asanyarray(n_spectra), np.asanyarray(date_passed), np.asanyarray(peak_mag)

data_table = Table([ZTF_id, ra, dec, z, classification, dist_to_galaxy, n_spectra, date_passed, peak_mag], names=("ZTF", "ra", "dec", "z", "class", "d_galaxy", "n_spec", "date_passed", "peak_app_mag"))

data_table.write("CLU_2019_12_18_backup.ascii", format='ascii') # saved new one 09/30

print ("Backup data has been completed! Seems like everything is running smooth... Now going to calculate the peak absolute magnitude")

############## QUERY HAS BEEN COMPLETED! ##############

# For each z value calculate the absolute_magnitude_peak
# For each z value calculate the absolute_magnitude_peak
for j in zip(data_table['z'], data_table['peak_app_mag']):
    z_start = j[0]
    pm_start= j[1]

    if z_start=="None" or pm_start=="None":
        print ("Sorry this redshift or peak magnitude is a None!")
        peak_abs_mag.append("None")
    else:
        print ("Okay looks like your values are free from NaN's")
        z = z_start.astype(float)
        pm = pm_start.astype(float)

        dist = redshift_to_distance(z).value
        print ("Here is the distance I calculated: %s"%dist)

        abs_mag = dist_mod_mag(pm, dist)
        print ("Here is the Abs_mag I calculated: %s"%abs_mag)
        peak_abs_mag.append(abs_mag)

print ("Here is the total number of absolute-magniutde peaks I found: %s"%len(peak_abs_mag))
print ("All looks good! Now finalizing the data-table for Census of the Local Universe!")
peak_abs_mag = np.asanyarray(peak_abs_mag)
#print (peak_abs_mag)

ZTF_id, ra, dec, z, classification = np.asanyarray(list(data_table['ZTF'])), np.asanyarray(list(data_table['ra'])), np.asanyarray(list(data_table['dec'])), np.asanyarray(list(data_table['z'])), np.asanyarray(list(data_table['class']))
dist_to_galaxy, n_spectra, date_passed, peak_mag = np.asanyarray(list(data_table["d_galaxy"])), np.asanyarray(list(data_table["n_spec"])), np.asanyarray(list(data_table["date_passed"])), np.asanyarray(list(data_table["peak_app_mag"]))

data_table_final = Table([ZTF_id, ra, dec, z, classification, dist_to_galaxy, n_spectra, date_passed, peak_mag, peak_abs_mag], names=("ZTF", "ra", "dec", "z", "class", "d_galaxy", "n_spec", "date_passed", "peak_app_mag", "peak_abs_mag"))

data_table_final.write("CLU_2019_12_18_original.ascii", format='ascii') # saved

print ("Congrats! You have queried all sources that have been stored on the CLU GROWTH Marshall.")
'''
