#!/usr/bin/python3

#
# The MIT License (MIT)
#
# Copyright (c) 2026 Michael J. Wouters
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

import argparse
import calendar
from datetime import datetime
from datetime import timezone
import os
import sqlite3
import sys
import time

# This is where cggttslib is installed
sys.path.append("/usr/local/lib/python3.8/site-packages")  # Ubuntu 20.04
sys.path.append("/usr/local/lib/python3.10/site-packages") # Ubuntu 22.04
sys.path.append("/usr/local/lib/python3.12/site-packages") # Ubuntu 24.04

try:
	import ottplib   as ottp
except ImportError:
	sys.exit('ERROR: Must install ottplib\n eg openttp/software/system/installsys.py -i ottplib')

VERSION = '0.2.0'
AUTHORS = 'Michael Wouters'

NPREVDAYS = 60

# ------------------------------------------
def ShowVersion():
	print (os.path.basename(sys.argv[0])+" "+VERSION)
	print ('Written by ' + AUTHORS)
	return

# ----------------------------------------------------------------------------------
home =os.environ['HOME']
root = home
configFile = os.path.join(root,'etc/butc.conf')
tmpDir  = os.path.join(root,'tmp')
nPrevDays = NPREVDAYS

if ottp.LibMajorVersion() >= 0 and ottp.LibMinorVersion() < 2: 
	sys.exit('Need ottplib minor version >= 2')

examples='TO DO'
parser = argparse.ArgumentParser(description='Releases data newly updated from Circular T for publication',
	formatter_class=argparse.RawDescriptionHelpFormatter,epilog=examples)
parser.add_argument('--config','-c',help='use an alternate configuration file',default=configFile)
parser.add_argument('--debug','-d',help='debug (to stderr)',action='store_true')
parser.add_argument('--version','-v',help='show version and exit',action='store_true')
parser.add_argument('--month','-m',help='release data for the given month 1..12 (current year assumed)')
parser.add_argument('--year','-y',help='release data for the given year (need to specify month too!)')

args = parser.parse_args()

if (args.version):
	ShowVersion()
	exit()

if (args.config):
	configFile = args.config
	if (not os.path.isfile(configFile)):
		ottp.ErrorExit(configFile + ' not found')

debug = args.debug
ottp.SetDebugging(debug)


cfg=ottp.Initialise(configFile,['main:gnss'])

if 'main:root' in cfg:
	root = cfg['main:root']
	
if ('database:file' in cfg):
	db = ottp.MakeAbsoluteFilePath(cfg['database:file'],root,os.path.join(home,'database'))


mjdToday  = ottp.MJD(time.time())
stopMJD   = mjdToday - 1 # previous day
startMJD  = stopMJD - nPrevDays

dt   = datetime.now(tz=timezone.utc) # get date in UTC
yyyy = dt.year

if args.year:
	if not args.month:
		ottp.ErrorExit('You need to specify the month as well (--month 1..12)')

if args.month: # manually specified a single month 
	if args.year:
		yyyy = int(args.year)
	mm = int(args.month)
	mmStart = mm
	yyyyStart = yyyy
	mmStop = mm
	yyyyStop = yyyy
	_,lastDayOfMonth = calendar.monthrange(yyyyStop,mmStop)
	startMJD = ottp.MJD(datetime(yyyyStart,mmStart,1,0,0,0,tzinfo=timezone.utc).timestamp())
	stopMJD =  ottp.MJD(datetime(yyyyStop,mmStop,lastDayOfMonth,0,0,0,tzinfo=timezone.utc).timestamp())

ottp.Debug(f'Processing for {startMJD} to {stopMJD}')

ottp.Debug(f'Connecting to the database {db}')
dbc = sqlite3.connect(db) # this opens a connection to the database
curs = dbc.cursor()  

# Identify the MJDs to be released
mjds=[]
for m in range(startMJD,stopMJD+1):
	r = curs.execute(f'SELECT release_utc FROM butcgnss WHERE mjd={m};')
	x = r.fetchone()
	if x:
		if x[0]==0:
			mjds.append(m)

if not mjds:
	print('Nothing to do')
	curs.close()
	dbc.close()
	sys.exit(0)
	
# Usually, there will be about 35
# so write 5 per line
print('Embargoed MJDs with UTC data that will be released for publication:\n')
for m in range(0,len(mjds)):
	print(mjds[m],end = ' ')
	if (m + 1) % 5 == 0:
		print() # end the line
if len(mjds) % 5:
	print('\n')
	
# Get confirmation
while True:
	ans = input("Release? (y/n): ").lower().strip()
	if ans in ['y', 'yes']:
		print('Releasing data ...')
		break
	elif ans in ['n', 'no']:
		print("Data release cancelled")
		curs.close()
		dbc.close()
		sys.exit(0)
	else:
		print("Please enter 'y' or 'n'.")
        
# Release the data
for m in mjds:
	cmd = f'UPDATE butcgnss SET release_utc = 1 WHERE mjd = {m}'
	curs.execute(cmd)

dbc.commit()  # and commit

# Release access to the DB so monthly report generation can run
curs.close()
dbc.close()

# Update reports
print('Updating the monthly reports ...')

# Push upload
print('Uploading to server ...')




