from lxml.html import parse, submit_form
from urllib2 import urlopen
from astropy.coordinates import SkyCoord
import astropy.units as u
from datetime import datetime
import argparse
import math
def mpcheck(ra,dec,searchrad=20,time=datetime.now()):
	coords = SkyCoord(ra=ra,dec=dec,unit=(u.degree,u.degree))
	ras = [str(int(coords.ra.hms[0])),str(int(coords.ra.hms[1])),'%.2f'%(coords.ra.hms[2])]
	print(ras)
	decs = [str(int(coords.dec.dms[0])),str(int(abs(coords.dec.dms[1]))),'%.2f'%(abs(coords.dec.dms[2]))]
	for i in range(len(ras)):
		if (float(ras[i])<10):
			ras[i] = '0' + ras[i]
			
	if (abs(float(decs[0]))<10):
		if (dec<0):
			decs[0] = '-0'+decs[0][1:]
		else :
			decs[0] = '0' +decs[0]
	for i in range(1,len(decs)):
		if (float(decs[i])<10):
			decs[i] = '0' + decs[i]
	fullra = ras[0] + ' ' + ras[1] + ' ' + ras[2]
	fulldec = decs[0] + ' ' + decs[1] + ' ' + decs[2]
	print ('Querying around %s %s at time %s'%(fullra,fulldec,str(time)))
	page = parse(urlopen("http://www.minorplanetcenter.net/cgi-bin/checkmp.cgi")).getroot()
	form = page.forms[1]
	form.fields['decl']= fulldec    
	form.fields['ra'] = fullra
	form.fields['type'] = 'p'
	form.fields['year'] = str(time.year)
	form.fields['month'] = str(time.month)
	form.fields['day'] = str(time.day + time.hour/24.)
	form.fields['radius'] = str(searchrad)
	#form.fields['radius'] = str(searchrad)
	olink = submit_form(form)
	#Form submitted
	opage = parse(olink).getroot()
	ochild = opage.getchildren()
	#obody = ochild[1].getchildren()
	opre = ochild[1].find('pre')
	ompc = opre.text_content() 
	#ompc = ompc.replace('\n')
	print(ompc)
	p = ompc.split('\n')
	p = p[4:-1]
	with open('mpcheck_ra_%s_dec_%s.reg'%(ra,dec),'w') as f:
		f.write('wcs\n')
		for i in p:
			coords = SkyCoord(ra=i[25:35],dec=i[36:47],unit=(u.hour,u.degree))
			f.write('point (%s,%s) # point=cross\n'%(coords.ra.degree,coords.dec.degree))
		f.close()
	f = open("mpcheck_objects_%s_%s"%(ra,dec),'w')
	f.write(ompc.encode('utf-8'))
	return(ompc)

if (__name__ == '__main__'):
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("ra",help="RA in degrees")
	parser.add_argument("dec",help="Dec in degrees")
	parser.add_argument("--searchrad",help="Search radius in arcmin", default=20)
	parser.add_argument("--time",help="Time at which query is wanted ex. 2019-02-27.32, default is now.")
	args = parser.parse_args()
	if args.time:
		t = args.time.split('-')
		th,td = math.modf(float(t[2])) 
		time = datetime(year=int(t[0]),month=int(t[1]),day=int(td),hour=int(round(th,2)*24))
	else:
		time = datetime.now()
	mpcheck(float(args.ra),float(args.dec),float(args.searchrad),time)