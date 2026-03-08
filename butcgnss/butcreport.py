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
import sqlite3
import sys

try:
	import ottplib   as ottp
except ImportError:
	sys.exit('ERROR: Must install ottplib\n eg openttp/software/system/installsys.py -i ottplib')

VERSION = '0.0.1'
AUTHORS = 'Michael Wouters'

MIN_UTCK_U = 5
MIN_UTC_U  = 6

# ------------------------------------------
def ShowVersion():
	print (os.path.basename(sys.argv[0])+" "+VERSION)
	print ('Written by ' + AUTHORS)
	return

# ---------------------------------------------
def WriteHeader(fout,header,mm,yyyy):
	
	if (os.path.isfile(header)):
		with open(header, 'r') as fin:
			txt = fin.read()
			fout.write('\n\n')
			fout.write(txt)
	else:
		ottp.Debug(f'Unable to open {header}')
	
	dt = datetime(yyyy,mm,1)
	fout.write('Report for {}\n'.format(dt.strftime('%B %Y')))
	
	# Date column is 10x
	# MJD  column is 5x
	# Each GNSS column is 9+5+9+5 = 28 x
	hdr1 = ' '*18
	hdr2 =  hdr4 = ' '*17 
	hdr3 = ' '*11 + '  MJD '
	
	for gi in range(0,len(gnss)):
		g=gnss[gi]
		hdr1 += '-'*12 + g + '-'*12
		
		hdr3 += '{:>9}{:>5}{:>9}{:>5}'.format('UTC('+lab+')','u','UTC','u')
		hdr4 += '{:>9}{:>14}'.format('- bUTC','- bUTC')
		if gi < len(gnss):
			hdr1 += ' '*2
			hdr3 += ' '
			hdr4 += ' '*6
	fout.write(hdr1+'\n')
	fout.write(hdr3+'\n')
	fout.write(hdr4+'\n')
	
	hdr = '-'*(16 + len(gnss)*28 + len(gnss))
	fout.write(hdr+'\n')
		
# ---------------------------------------------
def WriteFooter(fout,footer):
	if (os.path.isfile(footer)):
		with open(footer, 'r') as fin:
			txt = fin.read()
			fout.write('\n\n')
			fout.write(txt)
	else:
		ottp.Debug(f'Unable to open {footer}')
		
# ----------------------------------------------------------------------------------
home =os.environ['HOME']
root = home

configFile = os.path.join(root,'etc/butc.conf')
tmpDir  = os.path.join(root,'tmp')

footer = os.path.join(root,'etc/butcgnss_footer.txt')
header = os.path.join(root,'etc/butcgnss_header.txt')

if ottp.LibMajorVersion() >= 0 and ottp.LibMinorVersion() < 2: 
	sys.exit('Need ottplib minor version >= 2')

examples='TO DO'
parser = argparse.ArgumentParser(description='Updates the monthly report',
	formatter_class=argparse.RawDescriptionHelpFormatter,epilog=examples)
parser.add_argument('--config','-c',help='use an alternate configuration file',default=configFile)
parser.add_argument('--debug','-d',help='debug (to stderr)',action='store_true')
parser.add_argument('--version','-v',help='show version and exit',action='store_true')
parser.add_argument('--month','-m',help='create a report for the givem month 1..12 (current year assumed)')
parser.add_argument('--year','-y',help='create a report for the given year (need to specify month too!)')

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
gnss = [g.strip() for g in gnss] 

if 'main:root' in cfg:
	root = cfg['main:root']

if 'main:lab' in cfg:
	lab = cfg['main:lab']

if 'report:footer' in cfg:
	footer = ottp.MakeAbsoluteFilePath(cfg['report:footer'],root,footer)

if 'report:header' in cfg:
	header = ottp.MakeAbsoluteFilePath(cfg['report:header'],root,header)
	
if ('database:file' in cfg):
	db = ottp.MakeAbsoluteFilePath(cfg['database:file'],root,os.path.join(home,'database'))

dt   = datetime.now(tz=timezone.utc) # get date in UTC
mm   = dt.month
yyyy = dt.year
		
if args.year:
	if not args.month:
		ottp.ErrorExit('You need to specify the month as well (--month 1..12)')

if args.month: # manually specified 
	if args.year:
		yyyy = int(args.year)
	mm = int(args.month)
	
startMJD = ottp.MJD(datetime(yyyy,mm,1,0,0,0,tzinfo=timezone.utc).timestamp())
stopMJD =  startMJD + calendar.monthrange(yyyy,mm)[1] - 1

ottp.Debug(f'Processing for {startMJD} to {stopMJD}')

ottp.Debug(f'Connecting to the database {db}')
dbc = sqlite3.connect(db) # this opens a connection to the database
curs = dbc.cursor()   

missingData = '*' # FIXME configurable ?

minUncerts = {}
for g in gnss:
	gl = g.lower()
	minUncerts[g] = [MIN_UTCK_U,MIN_UTC_U]
	if gl+':utck uncertainty' in cfg:
		minUncerts[g][0] = float(cfg[gl+':utck uncertainty'])
	if gl+':utc uncertainty' in cfg:
		minUncerts[g][1] = float(cfg[gl+':utc uncertainty'])
		
for m in range(startMJD,stopMJD+1):
	
	ts = (m - 40587)*86400 # convert MJD to UNIX time
	dt = datetime.fromtimestamp(ts, tz=timezone.utc) # get date in UTC
	
	outputLine = f'{dt.year:4d}-{dt.month:02d}-{dt.day:02d} {m:5d} ' # 17 characters
	
	for gi in range(0,len(gnss)):
		g = gnss[gi]
		r = curs.execute(f'SELECT UTCk_{g},UTCk_{g}_u,UTC_{g},UTC_{g}_u from butcgnss where MJD={m};')
		x = r.fetchone()
		
		if x:
			minUTCkUncertainty = minUncerts[g][0]*10
			minUTCUncertainty  = minUncerts[g][1]
			
			if x[0]==None: # no data
				outputLine += f'{missingData:>9}{missingData:>5}{missingData:>9}{missingData:>5}'
			elif x[2]==None: # no UTC data
				reportedUTCkUncert = x[1]
				if reportedUTCkUncert  < minUTCkUncertainty:
					reportedUTCkUncert = minUTCkUncertainty
				outputLine += '{:>9}{:>5}{:>9}{:>5}'.format(round(x[0],1),math.ceil(reportedUTCkUncert),missingData,missingData)
			else: # Yay got it all
				reportedUTCkUncert = x[1]
				if reportedUTCkUncert  < minUTCkUncertainty:
					reportedUTCkUncert = minUTCkUncertainty
				reportedUTCUncert = x[3]
				if reportedUTCUncert  < minUTCUncertainty:
					reportedUTCUncert = minUTCUncertainty
				#ottp.Debug(f'{g}: UTC - bUTC {x[2]} +/- {x[3]}') # easiest place to put this
				outputLine += '{:>9}{:>5}{:>9}{:>5}'.format(round(x[0],1),math.ceil(reportedUTCkUncert),round(x[2],1),math.ceil(reportedUTCUncert))
		else:
			outputLine += f'{missingData:>9}{missingData:>5}{missingData:>9}{missingData:>5}'
		if gi < len(gnss):
				outputLine += ' '
				
	print(outputLine)

curs.close()
dbc.close()

currReport = None
if currReport:
	WriteFooter(currReport,footer)
	currReport.write('\n\n####################################################################\n')
	currReport.write(f'Generated by {os.path.basename(sys.argv[0])} v{VERSION}\n')
	if (mjdLastUTCupdate == mjd): # FIXME may be wrong
		currReport.write('UTC - bUTC_GNSS updated {}\n'.format(datetime.now(tz=timezone.utc).strftime('%Y-%m-%d')))
	currReport.write('####################################################################\n')
	currReport.close()		
	

		
