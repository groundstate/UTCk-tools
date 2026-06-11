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
import math
import os
from pathlib import Path
import sqlite3
import sys

# This is where cggttslib is installed
sys.path.append("/usr/local/lib/python3.8/site-packages")  # Ubuntu 20.04
sys.path.append("/usr/local/lib/python3.10/site-packages") # Ubuntu 22.04
sys.path.append("/usr/local/lib/python3.12/site-packages") # Ubuntu 24.04

try:
	import ottplib   as ottp
except ImportError:
	sys.exit('ERROR: Must install ottplib\n eg openttp/software/system/installsys.py -i ottplib')

VERSION = '0.1.0'
AUTHORS = 'Michael Wouters'

# ------------------------------------------
def ShowVersion():
	print (os.path.basename(sys.argv[0])+" "+VERSION)
	print ('Written by ' + AUTHORS)
	return
	
# ----------------------------------------------------------------------------------
home =os.environ['HOME']
root = home
lab = 'AUS'

configFile = os.path.join(root,'etc/butc.conf')

examples='EXAMPLES TO DO'
parser = argparse.ArgumentParser(description='Utility for eg querying the database',
	formatter_class=argparse.RawDescriptionHelpFormatter,epilog=examples)
parser.add_argument('mjd',nargs = '*',help='first MJD [last MJD] (if not given, the MJD of the previous day is used as the last MJD)')
parser.add_argument('--month','-m',help='specify month 1..12' )
parser.add_argument('--year','-y',help='specify year')
parser.add_argument('--gnss','-g',help='comma-separated list of GNSS to extract (BDS,GAL,GLO,BDS)')
parser.add_argument('--config','-c',help='use an alternate configuration file',default=configFile)
parser.add_argument('--debug','-d',help='debug (to stderr)',action='store_true')
parser.add_argument('--version','-v',help='show version and exit',action='store_true')

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

gnss = cfg['main:gnss'].split(',')
if args.gnss:
	gnss = args.gnss.split(',')
gnss = [g.strip() for g in gnss] 

if 'paths:root' in cfg:
	root = cfg['paths:root']
	# If the root is not absolute, prepend the user's home directory
	tmpPath = Path(root)
	if not tmpPath.is_absolute():
		root = os.path.join(home,root)

ottp.Debug(f'root path = {root}')
	
tmpDir  = os.path.join(root,'tmp')
if 'paths:tmp' in cfg:
	tmpDir = ottp.MakeAbsolutePath(cfg['paths:tmp'],root)
	
if 'main:lab' in cfg:
	lab = cfg['main:lab']
	
db = os.path.join(root,'butcgnss.db')
if ('database:file' in cfg):
	db = ottp.MakeAbsoluteFilePath(cfg['database:file'],root,os.path.join(home,'database'))

# If an MJD range is manually specified then
# processing is restricted to that range
if (args.mjd):
	
	if 1 == len(args.mjd):
		startMJD = int(args.mjd[0])
		stopMJD  = startMJD
	elif 2 == len(args.mjd):
		startMJD = int(args.mjd[0])
		stopMJD  = int(args.mjd[1])
		if (stopMJD < startMJD):
			ottp.ErrorExit('Stop MJD is before start MJD')
	else:
		ottp.ErrorExit('Too many MJDs')

if (args.month is not None) != (args.year is not None):
	ottp.ErrorExit('You need to specify --month and --year')
elif (args.month and args.year):
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

missingData = '*' # FIXME configurable ?

for m in range(startMJD,stopMJD+1):
		
	outputLine = f'{m:5d} ' # 17 characters
	
	for gi in range(0,len(gnss)):
		g = gnss[gi]
		r = curs.execute(f'SELECT UTCk_{g},UTCk_{g}_u,UTC_{g},UTC_{g}_u,release_utc from butcgnss where MJD={m};')
		x = r.fetchone()
		
		if x:
			if x[0]==None: # no data
				outputLine += f'{missingData:>}{missingData:>}{missingData:>}{missingData:>}'
			elif x[2] == None : # no UTC data
				outputLine += '{:>} {:>} {:>} {:>}'.format(x[0],x[1],missingData,missingData)
			else: # Yay got it all
				outputLine += '{:>} {:>} {:>} {:>}'.format(x[0],x[1],x[2],x[3])
		else:
			outputLine += f'{missingData:>}{missingData:>}{missingData:>}{missingData:>}'
		if gi < len(gnss):
				outputLine += ' '

		print(outputLine)
		
curs.close()
dbc.close()

	
	

		
