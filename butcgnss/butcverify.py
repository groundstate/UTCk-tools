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
import json
import numpy as np
import os
import requests
import sqlite3
import sys

try:
	import ottplib   as ottp
except ImportError:
	sys.exit('ERROR: Must install ottplib\n eg openttp/software/system/installsys.py -i ottplib')

VERSION = '0.0.1'
AUTHORS = 'Michael Wouters'

# ------------------------------------------
def ShowVersion():
	print (os.path.basename(sys.argv[0])+" "+VERSION)
	print ('Written by ' + AUTHORS)
	return

# ---------------------------------------------
def FetchBIPMBUTC(gnss,startMJD,stopMJD):
	
	ottp.Debug(f'Fetching Circular bUTC prediction data  for {gnss} from BIPM');
	#req = f'{httpRequest}scale=b_{gnss}&mjd1={startMJD}&mjd2={stopMJD}&outfile=json'
	#ottp.Debug(req)
	#try:
		#r = requests.get(req,verify=rootCert)
	#except:
		#ottp.Debug('http request failed')
		#return None,None,None
	
	#d = json.loads(r.text)
	#if not(d['errorcode'] == 0):
		#ottp.Debug(f"http request return error {d['errorcode']}")
		#return None,None,None
	
	#TEMPORARY CODE
	with open('/home/michael/database/butc_gps.json', 'r', encoding='utf-8') as file:
		d = json.load(file)
    
	#data = [[m,b,u] for m,b,u in zip(d['data'][0]['x'],d['data'][0]['y'],d['data'][0]['unc'])]

	return np.array(d['data'][0]['x']),np.array(d['data'][0]['y']),np.array(d['data'][0]['unc'],dtype=float)

def ReadLabBUTC(db,g,startMJD,stopMJD):
	mjd = []
	d   = []
	u   = []
	ottp.Debug(f'Connecting to the database {db}')
	dbc = sqlite3.connect(db) # this opens a connection to the database
	curs = dbc.cursor()   

	for mjd in range(startMJD,stopMJD+1):
		r = curs.execute(f'SELECT UTC_{g},UTC_{g}_u,release_utc from butcgnss where MJD={mjd};')
		x = r.fetchone()
		if x:
			if not (x[0] == None):
				mjd.append(m)
				d.append(x[0])
				u.append(x[1])
	curs.close()
	dbc.close()
	return np.array(m),np.array(d),np.array(u)
	 
# ----------------------------------------------------------------------------------
home =os.environ['HOME']
root = home

httpRequest = 'https://webtai.bipm.org/api/v1.0/get-data.html?'
rootCert = None # if you use an empty string, this will skip SSL verification which is a bad thing

configFile = os.path.join(root,'etc/butc.conf')
tmpDir  = os.path.join(root,'tmp')

if ottp.LibMajorVersion() >= 0 and ottp.LibMinorVersion() < 2: 
	sys.exit('Need ottplib minor version >= 2')

examples='TO DO'
parser = argparse.ArgumentParser(description='Provides information for checking a monthly report',
	formatter_class=argparse.RawDescriptionHelpFormatter,epilog=examples)
parser.add_argument('--config','-c',help='use an alternate configuration file',default=configFile)

parser.add_argument('--debug','-d',help='debug (to stderr)',action='store_true')
parser.add_argument('--display',help='display plots',action='store_true')
parser.add_argument('--version','-v',help='show version and exit',action='store_true')
parser.add_argument('--month','-m',help='create a report for the given month 1..12 (current year assumed)')
parser.add_argument('--year','-y',help='create a report for the given year (need to specify month too!)')

args = parser.parse_args()

if args.display:
	import matplotlib.pyplot as plt
else:
	import matplotlib as mplt # this (and the next line) stops warnings about being unable to connect to a display
	mplt.use('Agg') 
	import matplotlib.pyplot as plt
	
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

if ('database:file' in cfg):
	db = ottp.MakeAbsoluteFilePath(cfg['database:file'],root,os.path.join(home,'database'))


dt   = datetime.now(tz=timezone.utc) # get date in UTC
mm   = dt.month
yyyy = dt.year
dd   = dt.day

mmStart = mm - 2 # need to check the previous two months
yyyyStart = yyyy
if mmStart <= 0:
	mmStart += 12
	yyStart = yyyy - 1
	
mmStop = mm
yyyyStop = yyyy
lastDayOfMonth = dd -1 # for the current month, yesterday

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

ottp.Debug(f'{mmStart} {yyyyStart},{mmStop} {yyyyStop}')

startMJD = ottp.MJD(datetime(yyyyStart,mmStart,1,0,0,0,tzinfo=timezone.utc).timestamp())
stopMJD =  ottp.MJD(datetime(yyyyStop,mmStop,lastDayOfMonth,0,0,0,tzinfo=timezone.utc).timestamp())

ottp.Debug(f'Processing for {startMJD} to {stopMJD}')

if ('main:root certificate' in cfg):
	rootCert= cfg['main:root certificate']


ottp.Debug(f'Connecting to the database {db}')
dbc = sqlite3.connect(db) # this opens a connection to the database
curs = dbc.cursor()   

# Make a plot for each GNSS

fig, axs = plt.subplots(len(gnss),1,figsize=[8,12],squeeze=False) #unsqueeze so we always get an array of axes
dstr = f'{mmStop} {yyyyStop}'
title = 'bUTC_GNSS check for ' + dstr + '\n'
title += os.path.basename(sys.argv[0])+ ' v' + VERSION   + ' run ' + dt.strftime('%Y-%m-%d %H:%M:%S') + '\n'
fig.suptitle(title,ha='left',x=0.1)

plotIndex = 0

for g in gnss:
	m,d,u = FetchBIPMBUTC(g,startMJD,stopMJD)  # nb u is string!
	
	axs[plotIndex,0].plot(m,d,color='g')
	axs[plotIndex,0].scatter(m,d,s=9,marker='o',color='g')
	axs[plotIndex,0].plot(m,d+u,linestyle='dashed',color='g') # uncertainty
	axs[plotIndex,0].plot(m,d-u,linestyle='dashed',color='g') # uncertainty
	axs[plotIndex,0].set_ylabel('ns')
	axs[plotIndex,0].set_xlabel('MJD')
	axs[plotIndex,0].set_title(f'UTC - bUTC_{g}')
	axs[plotIndex,0].grid()
	
	m,d,u = ReadLabBUTC(db,g,startMJD,stopMJD)
	
	plotIndex += 1



if args.display:
	plt.show()

