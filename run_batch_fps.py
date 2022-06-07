from astropy.io import ascii
import os

t = ascii.read('/Users/viraj/ZTF/clu_subluminous_events_stats.csv')
inte = t[t['true_classification'] == '?']
inte = inte[inte['spec_exists']=='no']

inte = inte[51:]
print('Will submit request for %s sources'%(len(inte)))

for src in inte:
	
	print('python /Users/viraj/ztf_utils/submit_fps_request.py %s --fritz'%(src['Name']))
	os.system('python /Users/viraj/ztf_utils/submit_fps_request.py %s --fritz'%(src['Name']))
