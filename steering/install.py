#!/usr/bin/python3

import os
import pathlib
import shutil

home =os.environ['HOME'] 

# Make required directories

dirs = ['control','control/processed_steer','control/scheduled_steer','tmp','reports','etc']
for d in dirs:
	p = os.path.join(home,d)
	if not(os.path.exists(p)):
		print('Creating ' + p)
		os.mkdir(p)
	else:
		print(p + ' already exists')
		
# Make required files

files = ['control/ENABLE_STEERING','control/LAST_UTCR_DOWNLOAD','control/mask.dat']
for f in files:
	p = os.path.join(home,f)
	if not(os.path.exists(p)):
		print('Creating ' + p)
		pathlib.Path(p).touch()
	else:
		print(p + ' already exists')
		
# Copy config file across
cfgdir = os.path.join(home,'etc')
shutil.copyfile('utcsteersched.conf',os.path.join(cfgdir,'utcsteersched.conf'))
