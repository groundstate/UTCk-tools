#!/usr/bin/python3

#
# The MIT License (MIT)
#
# Copyright (c) 2025 Michael J. Wouters
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
import re
import requests
import shutil
import sys
import time

# This is where cggttslib is installed
sys.path.append("/usr/local/lib/python3.8/site-packages")  # Ubuntu 20.04
sys.path.append("/usr/local/lib/python3.10/site-packages") # Ubuntu 22.04
sys.path.append("/usr/local/lib/python3.12/site-packages") # Ubuntu 24.04

try: 
	import cggttslib as cggtts
except ImportError:
	sys.exit('ERROR: Must install cggttslib\n eg openttp/software/system/installsys.py -i cggttslib')

from cggttslib import CGGTTS

try:
	import ottplib   as ottp
except ImportError:
	sys.exit('ERROR: Must install ottplib\n eg openttp/software/system/installsys.py -i ottplib')

try:
	import rinexlib   as rinex
except ImportError:
	sys.exit('ERROR: Must install rinexlib\n eg openttp/software/system/installsys.py -i rinexlib')
	
VERSION = '0.4.0'
AUTHORS = 'Michael Wouters'

SQRT2 = math.sqrt(2)

REFSYS_AVG_WINDOW = 2 # window size for averaging


# Uncertainty of GNSS provider's link calibration 
# Values from Defraigne et al 2023 Metrologia 60 065010,  DOI : 10.1088/1681-7575/ad0562
U_CAL_GNSS = {'BDS': 2.4, 'GAL': 2.4 ,'GLO': 3.8, 'GPS': 2.7}

# Uncertainty arising from choice of UTC model
# Values from Defraigne et al 2023 Metrologia 60 065010,  DOI : 10.1088/1681-7575/ad0562
U_NAVMSG_GNSS = {'BDS': 0.2, 'GAL': 0.1 ,'GLO':1.2, 'GPS': 1.3}

CIRT_UA = 0.2            # uA as reported in Circular T

CLOCK_5071_STD = 1 # standard model
CLOCK_5071_HPT = 2 # high performance model 
CLOCK_MASER    = 3

# ------------------------------------------
def ShowVersion():
	print (os.path.basename(sys.argv[0])+" "+VERSION)
	print ('Written by ' + AUTHORS)
	return

#----------------------------------------------
def CheckConfig(cfg,req):
	ok = True
	for r in req:
		if not(r in cfg):
			print(f'{r} is not set')
			ok = False
	return ok

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

# ------------------------------------------
def FindCGTTSFile(gnss,mjd):
	# Find the CGGTTS file
	fname = None
	gToken = gnss.lower()
	path = os.path.join(root,cggttsDir)
	if (gToken+':cggtts path') in cfg:
		path = ottp.MakeAbsolutePath(cfg[gToken+':cggtts path'],root)
	prefix = ''
	if (gToken+':cggtts prefix') in  cfg:
		prefix = cfg[gToken+':cggtts prefix']
	if prefix:
		ext = ''
	else:
		ext = 'cctf'
	if (gToken+':cggtts extension') in  cfg:
		ext = cfg[gToken+':cggtts extension']
	
	fname = cggtts.FindFile(path,prefix,ext,mjd)
	if (not os.path.isfile(fname)): 
		ottp.Debug(f"Couldn't find primary CGGTTS file path={path} prefix={prefix} extension={ext}")
		# Look for alternate
		if (gToken + ':alternate cggtts path') in cfg:
			prefix = ''
			if (gToken+':alternate cggtts prefix') in  cfg:
				prefix = cfg[gToken+':alternate cggtts prefix']
			ext = 'cctf'
			if (gToken+':alternate cggtts extension') in  cfg:
				ext = cfg[gToken+':alternate cggtts extension']
			fname = cggtts.FindFile(path,prefix,ext,mjd)
			if (not os.path.isfile(fname)): 
				ottp.Debug(f"Couldn't find alternate CGGTTS file path={path} prefix={prefix} extension={ext} - must be the End of Days")
	return fname
	
# Returns time system corrections as a dictionary, with the GNSS name as the key
# and the model parameters in the order of the RINEX file
# navFile is presumed to be openable
## TODO earlier RINEX ??
# ------------------------------------------
def GetTimeSysCorr(navFile):
	
	tsc = {'BDS': None, 'GAL': None ,'GLO': None, 'GPS': None}
	fnav = open(navFile,'r')
	for l in fnav:
		# Fields are offset a0, rate a1, tow, wn 
		# covers version >= 2.12 
		if l[0:4] == 'BDUT': # covers version 2.12 ->
			tsc['BDS'] = [float(l[5:22].replace('D','E')), float(l[22:38].replace('D','E')), int(l[38:45]),int(l[45:50])]
		elif l[0:4] == 'GAUT':
			tsc['GAL'] = [float(l[5:22].replace('D','E')), float(l[22:38].replace('D','E')), int(l[38:45]),int(l[45:50])]	
		elif  l[0:4] == 'GLUT':
			tsc['GLO'] = [float(l[5:22].replace('D','E')), float(l[22:38].replace('D','E')), int(l[38:45]),int(l[45:50])]	
		elif  l[0:4] == 'GPUT': 
			tsc['GPS'] = [float(l[5:22].replace('D','E')), float(l[22:38].replace('D','E')), int(l[38:45]),int(l[45:50])]
		elif  l[60:73] == 'END OF HEADER': # CHECKED
			break
	fnav.close()
	return tsc

# ------------------------------------------
def DeltaGNSSUTC(gnss,Wn,Dn,leapSecs,tsCorr):
	deltaUTC = 0
	if (gnss == 'GPS') or (gnss == 'GAL'):
		# From the GPS ICD
		# t_UTC = t_E - delta_UTC where t_E is 'effective' GPS time 
		# delta_UTC = dt_LS + A0 + A1*(t_E - t0t + 604800*(WN - Wn_t) )
		# We want the value at UTC0 for the day which means
		# t_E = 86400*gpsDn + leapSecs
		# GAL is the same
		deltaUTC = 1.0E9*(tsCorr[0] + tsCorr[1]*(86400*Dn + leapSecs - tsCorr[2] + 604800*(Wn - tsCorr[3]))) # in ns
	elif (gnss == 'GLO'):
		deltaUTC = 1.0E9*tsCorr[0]
	ottp.Debug(f'DeltaGNSSUTC: {gnss} {deltaUTC}')
	return deltaUTC

# ---------------------------------------------
def GetRefsys(cgBef,cgAft,winSize):

	wSz = 1800*winSize # winSize is in hours
	refsys0 = None
	uRefsys0 = None
	nTracks = 0
	
	tBef = tAft = 0
	nBef = nAft = 0

	# Seems messy compared with a concatenating tha data into a single array
	# it does allow the two CGGTTS files to have different formats
	
	if cgBef.fileName:
		tStart = 86400 - wSz
		ntrks = len(cgBef.tracks)
		for t in range(0,ntrks):
			if cgBef.tracks[t][cgBef.STTIME] >= tStart: # don't worry about the 390 s offset
				tBef = t 
				nBef = ntrks - tBef
				break
				
	if cgAft.fileName:
		tStop = wSz
		ntrks = len(cgBef.tracks)
		for t in range(0,ntrks):
			if cgAft.tracks[t][cgBef.STTIME] >= tStop:
				tAft = t - 1
				nAft = tAft + 1
				break
		
	if not(nBef and nAft):
		return refsys0,uRefsys0,nTracks
	
	refsys0 = 0
	
	if nBef > 0:
		ntrks = len(cgBef.tracks)
		for t in range(tBef,ntrks):
			#print(cgBef.tracks[t][cgBef.STTIME],cgBef.tracks[t][cgBef.REFSYS])
			refsys0 = refsys0 + cgBef.tracks[t][cgBef.REFSYS]
	if nAft > 0:
		ntrks = len(cgAft.tracks)
		for t in range(0,tAft+1):
			#print(cgAft.tracks[t][cgAft.STTIME],cgAft.tracks[t][cgAft.REFSYS])
			refsys0 = refsys0 + cgAft.tracks[t][cgAft.REFSYS]
	
	nTracks = nBef + nAft
	if nTracks == 1: # wayy... too little
		return refsys0,uRefsys0,nTracks
		
	refsys0 = refsys0/nTracks 
	
	urefsys0 = 0
	if nBef > 0:
		ntrks = len(cgBef.tracks)
		for t in range(tBef,ntrks):
			urefsys0 += (refsys0 - cgBef.tracks[t][cgBef.REFSYS])**2
	if nAft > 0:
		ntrks = len(cgAft.tracks)
		for t in range(0,tAft+1):
			urefsys0 += (refsys0  - cgAft.tracks[t][cgAft.REFSYS])**2
	
	return refsys0,math.sqrt(urefsys0/(nTracks-1)),nTracks  # returns std dev as uncertainty for refsys0

# Returns Circular T in a list
# ---------------------------------------------
def GetCircularT(lab,startMJD,stopMJD):
	# Code for testing - be kind to the BIPM web server
	data = []
	fin = open(os.path.join(root,'report/cirt.txt'),'r')
	firstMJD = lastMJD = None
	for l in fin:
		if l[0] == '#':
			continue
		vals = l.strip().split()
		if (len(vals) == 3):
			data.append([int(vals[0]),float(vals[1]),float(vals[2]),
				uACircularT,math.sqrt(float(vals[2])**2 - uACircularT**2 )]) #  TODO when the web API reports (uA,uB)
			if not firstMJD:
				firstMJD = int(vals[0])
			lastMJD = int(vals[0])
	fin.close()
	return data,firstMJD,lastMJD

# ---------------------------------------------
def __GetCircularT(lab,startMJD,stopMJD):
	ottp.Debug('Fetching Circular T');
	ottp.Debug(f'{httpRequest}scale=utc&lab={lab}&mjd1={startMJD}&mjd2={stopMJD}&outfile=txt');
	try:
		r = requests.get(f'{httpRequest}scale=utc&lab={lab}&mjd1={startMJD}&mjd2={stopMJD}&outfile=txt')
	except:
		return None,None,None
	
	lines = r.text.split('\r\n')
	data = []
	firstMJD = lastMJD = None
	for l in lines:
		l = l.strip()
		if not(l):
			continue
		if l[0] == '#':
			continue
		vals = l.strip().split()
		if (len(vals) == 3):
			data.append([int(vals[0]),float(vals[1]),float(vals[2]),
				uACircularT,math.sqrt(float(vals[2])**2 - uACircularT**2 )]) #  TODO when the web API reports (uA,uB)
			if not firstMJD:
				firstMJD = int(vals[0])
			lastMJD = int(vals[0])
	return data,firstMJD,lastMJD
	
# ---------------------------------------------
# This uses the list representation!
# 0->MJD, 1->UTC-UTCk, 2->u, 3->uA, 4->uB
def GetNearestCirtU(cirt,mjd):
	for i in range(0,len(cirt)-1):
		if (cirt[i][0] >= mjd and mjd <= cirt[i+1][0]):
			ottp.Debug(f'{cirt[i][0]} {cirt[i+1][0]} {mjd} {cirt[i][2]} ')
			return cirt[i][3], cirt[i][4]
	ottp.Debug(f'{cirt[-1][0]} {mjd} failed -> {cirt[i][3]}  {cirt[i][4]}')
	return cirt[-1][3],cirt[-1][4]

# ---------------------------------------------
# Linear interpolation  of Circular T
# Returns the interpolated value and its uncertainty
# This uses the dict representation!

def InterpolateUTC(mjd,mjd0,mjd1,cirt):
	utc0 = cirt[mjd0][1]
	utc1 = cirt[mjd1][1]	
	dmjd = mjd - mjd0
	
	# The uncertainty is calculated by
	# using the recommended uncertainty for the mean frequency, as given in CCTF WGMRA Guideline 4 
	uFreq = cirt[mjd0][3] * SQRT2/ 5.0  # fractional, ns/day
	minD = dmjd # use the minimum distance 
	if mjd1-mjd < minD:
		minD = mjd1-mjd
	utcm = utc0 + dmjd*(utc1-utc0)/5 # estimate of UTC - UTC(k) by linear interpolation
	uc   = math.sqrt((uFreq*minD)**2 ) 
	ottp.Debug(f'InterpolateUTC {mjd} {utcm} {uc}, [{utc0} {utc1}]')
	return utcm,uc


# -------------------------------------------
def ClockStability(clockModel,tauDays):
	if clockModel == CLOCK_5071_STD:
		return 2.0*math.sqrt(tauDays) # in ns 
	elif clockModel == CLOCK_MASER:
		return 0.0 # FIXME
	else:
		sys.exit('Requested clock model is unsupported')

	
# ---------------------------------------------
# Main 
# ---------------------------------------------

home =os.environ['HOME']
root = home

configFile = os.path.join(root,'etc/butcgnss.conf')
rnxDir = os.path.join(root,'rinex')
tmpDir  = os.path.join(root,'tmp')
reportDir = os.path.join(root,'reports')
cggttsDir = os.path.join(root,'cggtts')
footer = os.path.join(root,'etc/butcgnss_footer.txt')
header = os.path.join(root,'etc/butcgnss_header.txt')
lab = 'AUS'
httpRequest = 'https://webtai.bipm.org/api/v1.0/get-data.html?'

winSize = REFSYS_AVG_WINDOW

minUTCkUncertainty = 5
minUTCUncertainty = 6
clockModel = CLOCK_5071_STD
uACircularT = CIRT_UA

if ottp.LibMajorVersion() >= 0 and ottp.LibMinorVersion() < 2: 
	sys.exit('Need ottplib minor version >= 2')

examples='TO DO'
parser = argparse.ArgumentParser(description='Generate UTC(k) - bUTC_GNSS and UTC - bUTC_GNSS ',
	formatter_class=argparse.RawDescriptionHelpFormatter,epilog=examples)
parser.add_argument('--config','-c',help='use an alternate configuration file',default=configFile)
parser.add_argument('--utc','-u',help='generate UTC-bUTC_gnss differences for the previous month',action='store_true')
parser.add_argument('--month','-m',help='generate differences for the givem month 1..12 (current year assumed)')
parser.add_argument('--year','-y',help='generate differences for the given year (need to specify month too!)')
parser.add_argument('--debug','-d',help='debug (to stderr)',action='store_true')
parser.add_argument('--force','-f',help='force overwriting of an existing reort',action='store_true')
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
cggtts.SetWarnings(debug)

cfg=ottp.Initialise(configFile,['main:gnss'])

gnss = cfg['main:gnss'].split(',')
gnss = [g.strip() for g in gnss] 

if 'main:root' in cfg:
	root = cfg['main:root']

if 'main:lab' in cfg:
	lab = cfg['main:lab']
	
if 'main:report path' in cfg:
	reportDir = ottp.MakeAbsolutePath(cfg['main:report path'],root)
	
if 'main:footer' in cfg:
	footer = ottp.MakeAbsoluteFilePath(cfg['main:footer'],root,footer)

if 'main:header' in cfg:
	header = ottp.MakeAbsoluteFilePath(cfg['main:header'],root,header)
	
if 'main:utck uncertainty' in cfg:
	minUTCkUncertainty = float(cfg['main:utck uncertainty'])
	
if 'main:utc uncertainty' in cfg:
	minUTCUncertainty = float(cfg['main:utc uncertainty'])

if 'main:clock' in cfg:
	clk = cfg['main:clock'].lower()
	if clk == 'hp5071 std':
		clockModel = CLOCK_5071_STD
	elif clk == 'hp5071 hpt':
		clockModel = CLOCK_5071_HPT
	elif clk == 'maser':
		clockModel = CLOCK_MASER
	else:
		ottp.ErrorExit(f"Unknown clock: {cfg['main:clock']}")
	
if 'main:circular t ua' in cfg:
	uACircularT = float(cfg['main:circular t ua'])

if 'rinex:path' in cfg:
	rnxDir = ottp.MakeAbsolutePath(cfg['rinex:path'],root)
	# check other necessary things have been defined
	if not(CheckConfig(cfg,['rinex:station name','rinex:version'])):
		ottp.ErrorExit('Missing entries in configuration file')
	staName = cfg['rinex:station name']
	rnxVersion = int(cfg['rinex:version'])

dailyUpdate = False
UTCupdate   = False

if args.year:
	if not args.month:
		ottp.ErrorExit('You need to specify the month as well (--month 1..12)')
	
if args.utc:
	
	if args.month: # manually specified 
		if args.year:
			yyyy = int(args.year)
		else:
			dt = datetime.now(tz=timezone.utc) # get date in UTC
			yyyy = dt.year
		mm   = int(args.month)
		startMJD = ottp.MJD(datetime(yyyy,mm,1,0,0,0,tzinfo=timezone.utc).timestamp())
		stopMJD =  startMJD + calendar.monthrange(yyyy,mm)[1] - 1
	else:
		
		# Processing for the previous two months so we need to work out what that is
		dt = datetime.now(tz=timezone.utc) # get date in UTC
		yyyy = dt.year
		mm   = dt.month
		dd   = dt.day
		
		if (dd < 12): # won't be available yet
			sys.exit('Too early for Circular T') # we're done
		if (mm == 1):
			mm1 = 12
			mm2 = 11
			yyyy1 = yyyy2 = yyyy - 1
		elif (mm == 2):
			mm1 = 1
			yyyy1 = yyyy
			mm2 = 12
			yyyy2 = yyyy-1
		else:
			mm1 = mm - 1 
			mm2 = mm - 2
			yyyy1 = yyyy2 = yyyy
	
		startMJD = ottp.MJD(datetime(yyyy2,mm2,1,0,0,0,tzinfo=timezone.utc).timestamp())
		stopMJD  = startMJD + calendar.monthrange(yyyy2,mm2)[1] + calendar.monthrange(yyyy1,mm1)[1] - 1
	
	startMJD -= 1
	# TESTED January and February 
	# Check if Circular T is available 
	cirtAsList,firstMJD,lastMJD = GetCircularT(lab,startMJD - 7,stopMJD + 7) 
	if not cirtAsList:
		sys.exit("Couldn't get Circular T")
	# In UTC update mode it's convenient to have Circular T as a dictionary
	cirt = {}
	for c in cirtAsList:
		cirt[c[0]] = [c[0],c[1],c[2],c[3],c[4]]
	
	# Check that the required MJDs are there
	if not(firstMJD <= startMJD and stopMJD - 5  <= lastMJD ): # probably won't get the last day of the month
		sys.exit(f"Not all required MJDs are available in Circular T (wanted [{startMJD},{stopMJD} - 5], got [{firstMJD},{lastMJD}])")
	
	UTCupdate = True
	
else:
	ts = time.time() - 86400 # previous day
	# Set the default processing range
	stopMJD = ottp.MJD(ts)  
	

	dt = datetime.fromtimestamp(ts, tz=timezone.utc) # get date in UTC
	ts   = ts - 86400*(dt.day-1) # FIXME check month rollover
	startMJD = ottp.MJD(ts)-1 # need to get the previous day as well
	mjdToday = startMJD + 1
	
	dailyUpdate = True 
	ottp.Debug(f'Default processing range {startMJD} - {stopMJD}')
	
	# Need to get Circular T so that we know the current uncertainty uB
	# Unfortunately the Web API returns the total uncertainty.
	# In most cases, and for us, uB dominates so for the moment
	# we'll use the total uncertainty for uB. Typically, we won't have a matching uB
	# for an MJD, so we'll just be using the most recent anyway (noting that this will be fixed
	# later when the UTC update is done)

	# Initial guess on range of Circular T to ask for is
	# [startMJD -7, stopMJD + 7] say
	
	cirtStartMJD = startMJD - 7
	cirtStopMJD  = stopMJD  + 7
	
	# The most recent value could be 6 or 7 weeks ago (two weeks into month, previous month not yet available)
	# It doesn't hurt to get a bit more data than we need so assume the most recent value is 50 days old.
	
	if (mjdToday - 50 < cirtStartMJD):
		cirtStartMJD = mjdToday - 50

	cirtAsList,firstMJD,lastMJD = GetCircularT(lab,cirtStartMJD,cirtStopMJD)
	if not cirtAsList:
		sys.exit("Couldn't get Circular T")
	
	# Range checking is fraught
	
	
currReport = None # initially, the report file is closed or doesn't yet exist
mjdLastUTCupdate = None # for tracking if we updated UTC - buTC_GNSS for the whole month
prevCGGTTSFile = {}
prevCGGTTSFile['BDS'] = CGGTTS(None,None)
prevCGGTTSFile['GAL'] = CGGTTS(None,None)
prevCGGTTSFile['GLO'] = CGGTTS(None,None)
prevCGGTTSFile['GPS'] = CGGTTS(None,None)
newData = {}

mjd = startMJD - 1

while mjd < stopMJD :
	mjd += 1 # doing this here means not having to increment multiple times
	reportData = {} # all GNSS, for this MJD
	
	ts = (mjd - 40587)*86400 # convert MJD to UNIX time
	dt = datetime.fromtimestamp(ts, tz=timezone.utc) # get date in UTC
	yyyy = dt.year
	mm   = dt.month
	dd   = dt.day
	doy = int(dt.strftime('%j'))
	gpsWn,gpsDn = ottp.MJDtoGPSWeekDay(mjd)
	ottp.Debug(f'\nProcessing {mjd}: {yyyy}-{mm}-{dd}, WN={gpsWn} DN={gpsDn}')
	
	if (dt.day == 1) and (dailyUpdate or UTCupdate): # beginning of month
		# Close any open file
		if currReport:
			WriteFooter(currReport,footer)
			currReport.write('\n\n####################################################################\n')
			currReport.write(f'Generated by {os.path.basename(sys.argv[0])} v{VERSION}\n')
			if (mjdLastUTCupdate == mjd -1):
				currReport.write('UTC - bUTC_GNSS updated {}\n'.format(datetime.now(tz=timezone.utc).strftime('%Y-%m-%d')))
			currReport.write('####################################################################\n')
			currReport.close()
			
		reportName = os.path.join(reportDir,f'butc{mm:02d}{yyyy:4d}.txt')
		
		if UTCupdate:
			# Check whether the file needs processing
			if os.path.isfile(reportName):
				fin = open(reportName,'r')
				with fin:
					txt = fin.read()
				if 'UTC - bUTC_GNSS updated' in txt:
					if not args.force:
						ottp.Debug(f'{reportName} is up to date ... skipping')
						# Skip ahead
						# Can do this by changing startMJD to the end of the month
						startMJD = startMJD + (calendar.monthrange(yyyy,mm)[1] - 1)  + 1
						mjd = startMJD - 1 # because we will imediately increment it
						continue
					else:
						ottp.Debug(f'{reportName} exists and is up to date ... overwriting')
		
		ottp.Debug(f'Creating report {reportName}')
		try:
			currReport = open(reportName,"w")
		except:
			sys.exit(f"Couldn't open {reportName}")
		
		WriteHeader(currReport,header,mm,yyyy)
	
	# Find the navigation file
	navFile,zExt = rinex.FindNavigationFile(rnxDir,staName,yyyy,doy,rnxVersion,False) # don't exit if not found

	timeSysCorr = {}
	if navFile:
		ottp.Debug(f'Found {navFile}, compression = {zExt}')
		src = navFile + zExt
		srcBase = os.path.basename(src)
		dst = os.path.join(tmpDir,srcBase)
		# We may not own the file so we need to make a local copy
		shutil.copy(src,dst)
		navFile,zAlgorithm = rinex.Decompress(dst)
		timeSysCorr = GetTimeSysCorr(navFile)
		leapSecs = rinex.GetLeapSeconds(navFile,rnxVersion) # read twice so a bit inefficient but anyway
		os.unlink(navFile)
		if debug:
			print(timeSysCorr)
	else:
		ottp.Debug("Couldn't find a navigation file for {yyyy} {doy}")
		
	for g in gnss:
		ottp.Debug(f'Processing {g}')
		reportData[g] = [None,None,None,None] # Fields are UTC(k)-bUTC_GNSS, u, UTC - bUTC_GNSS, u]
		
		fName = FindCGTTSFile(g,mjd)
		
		if fName:
			ottp.Debug(f'Reading {fName}')
			cgf = CGGTTS(fName,mjd)
			cgf.Read()
		else:
			ottp.Debug('No CGGTTS file found')
			cgf = CGGTTS(None,mjd)
		
		# The first one is for the day before the start MJD so no more to do
		if (mjd==startMJD):
			prevCGGTTSFile[g] = cgf # save it 
			continue
		
		# Calculate the average REFSYS, if we can
		refsys0 = None
		if prevCGGTTSFile[g].fileName or cgf.fileName:
			refsys0,uRefsys0,nTracks = GetRefsys(prevCGGTTSFile[g],cgf,winSize)
			#print(refsys0,uRefsys0,nTracks)
		
		if (refsys0 == None):
			ottp.Debug('Insufficent CGGTTS data')
			prevCGGTTSFile[g] = cgf
			continue
	
		if not navFile: # check this before we move on
			prevCGGTTSFile[g] = cgf
			continue
			
		tsCorr = timeSysCorr[g]
		if not tsCorr:
			ottp.Debug(f'time sys corrections unavailable for {g}')
			prevCGGTTSFile[g] = cgf
			continue
			
		
		uAcirt,uBcirt = GetNearestCirtU(cirtAsList,mjd)
		deltaUTC = DeltaGNSSUTC(g,gpsWn,gpsDn,leapSecs,tsCorr)
			
		reportData[g][0] = refsys0 + deltaUTC
		reportData[g][1] = math.sqrt(uRefsys0**2 + uBcirt**2 + U_CAL_GNSS[g]**2 +  U_NAVMSG_GNSS[g]**2) # only uB is relevant here
		ottp.Debug('UTCk - bUTC {} {:g} {:g} {:g} {:g} {:g}'.format(g,mjd, reportData[g][0],reportData[g][1],refsys0,uRefsys0))
		if UTCupdate:
			mjdLastDigit = int(str(mjd)[-1])
			if mjdLastDigit < 4:
				mjd0 = mjd - mjdLastDigit - 1
				mjd1 = mjd - mjdLastDigit + 4
			else:
				mjd0 = mjd - mjdLastDigit + 4
				mjd1 = mjd - mjdLastDigit + 9
			if not(mjd0 in cirt) and not(mjd1 in cirt):
				ottp.Debug(f'{mjd0} and {mjd1} not in CirT')
				prevCGGTTSFile[g] = cgf
				continue
				
			if not(mjd1 in cirt): # special case: can't interpolate but can do mjd ==  mjd0 and get one more day
				ottp.Debug(f'{mjd1} not in CirT')
				if (mjd == mjd0):
					reportData[g][2] = cirt[mjd][1] + reportData[g][0]
					reportData[g][3] = math.sqrt(uRefsys0**2 +  uAcirt**2 + uBcirt**2 + U_CAL_GNSS[g]**2 +  U_NAVMSG_GNSS[g]**2) # no contribution from instability
					mjdLastUTCupdate = mjd
				prevCGGTTSFile[g] = cgf # this is correct - we only have mjd0, so we're done after tewsing for mjd==mjd0 
				continue
			
			# Don't interpolate on CircularT days - smaller uncertainty
			# FIXE this is unnecessary
			#if (mjd == mjd0 or mjd==mjd1):
				#reportData[g][2] = cirt[mjd][1] + reportData[g][0]
				#reportData[g][3] = math.sqrt(uRefsys0**2 +  uAcirt**2 + uBcirt**2 + U_CAL_GNSS[g]**2 +  U_NAVMSG_GNSS[g]**2) # no contribution from instability
				#mjdLastUTCupdate = mjd
				#prevCGGTTSFile[g] = cgf
				#continue
				
			utc0 = cirt[mjd0][1]
			utc1 = cirt[mjd1][1]	
			
			utcDiff,utcInterpUncert = InterpolateUTC(mjd,mjd0,mjd1,cirt)
			reportData[g][2] = utcDiff + reportData[g][0]
			# Uncertainty sources:
			# time transfer noise == uRefsys0
			# UTC interpolation uncertainty == utcUncert (which includes the link calibration uncertainty)
			# UTC(k) instability 
			# GNSS provider link's calibration uncertainty
			# UTC prediction 
			UTCkInstability = ClockStability(clockModel,mjd - mjd0)
			if  ClockStability(clockModel,mjd - mjd0) <  UTCkInstability:
				UTCkInstability = ClockStability(clockModel,mjd1- mjd)
			
			reportData[g][3] = math.sqrt(uRefsys0**2 + cirt[mjd0][3]**2 + cirt[mjd0][4]**2 + utcInterpUncert**2 +
					UTCkInstability**2+ U_CAL_GNSS[g]**2 +  U_NAVMSG_GNSS[g]**2) 
			
			mjdLastUTCupdate = mjd # keep track of the last MJD we update so that we can test if a month is complete 
			
		prevCGGTTSFile[g] = cgf # save it for the next time
	
	if (mjd == startMJD):
		continue
	
	# If we can't generate the data, mark it as missing'
	outputLine = f'{yyyy:4d}-{mm:02d}-{dd:02d} {mjd:5d} ' # 17 characters
	missingData = '*'

	for gi in range(0,len(gnss)):
		g = gnss[gi]
		if reportData[g][0]==None:
			outputLine += f'{missingData:>9}{missingData:>5}{missingData:>9}{missingData:>5}'
		elif reportData[g][2]==None:
			reportedUTCkUncert = reportData[g][1]
			if reportedUTCkUncert  < minUTCkUncertainty:
				reportedUTCkUncert = minUTCkUncertainty
			outputLine += '{:>9}{:>5}{:>9}{:>5}'.format(round(reportData[g][0],1),math.ceil(reportedUTCkUncert),missingData,missingData)
		else:
			reportedUTCkUncert = reportData[g][1]
			if reportedUTCkUncert  < minUTCkUncertainty:
				reportedUTCkUncert = minUTCkUncertainty
			reportedUTCUncert = reportData[g][3]
			if reportedUTCUncert  < minUTCUncertainty:
				reportedUTCUncert = minUTCUncertainty
			ottp.Debug(f'{g}: UTC - bUTC {reportData[g][2]} {reportData[g][3]}') # easiest place to put this
			outputLine += '{:>9}{:>5}{:>9}{:>5}'.format(round(reportData[g][0],1),math.ceil(reportedUTCkUncert),round(reportData[g][2],1),math.ceil(reportedUTCUncert))
		if gi < len(gnss):
			outputLine += ' '
	if (dailyUpdate or UTCupdate):
		currReport.write(outputLine+'\n')
	newData[mjd] = outputLine + '\n'

# Finish off any currently open file	
if currReport:
	WriteFooter(currReport,footer)
	currReport.write('\n\n####################################################################\n')
	currReport.write(f'Generated by {os.path.basename(sys.argv[0])} v{VERSION}\n')
	if (mjdLastUTCupdate == mjd): # FIXME may be wrong
		currReport.write('UTC - bUTC_GNSS updated {}\n'.format(datetime.now(tz=timezone.utc).strftime('%Y-%m-%d')))
	currReport.write('####################################################################\n')
	currReport.close()		

	

	
