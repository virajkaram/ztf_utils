import os
from astropy.io import ascii

sourcenames = ascii.read('lrn-ilrt-ZTFI.csv')

for source in sourcenames:
	print('Currently doing %s'%(source['Name']))
	os.system('python submit_fps_request.py %s --fritz'%(source['Name']))
	#print('python submit_fps_request.py %s --fritz'%(source['Name']))
