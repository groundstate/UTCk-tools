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
import sqlite3
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
	
VERSION = '0.9.1'
AUTHORS = 'Michael Wouters'

SQRT2 = math.sqrt(2)

REFSYS_AVG_WINDOW = 2 # window size for averaging

NPREVDAYS = 60 # number of days to regenerate prior to the current day (when not running with explicit MJD range)

# Uncertainty of GNSS provider's link calibration 
# Values from Defraigne et al 2023 Metrologia 60 065010,  DOI : 10.1088/1681-7575/ad0562
U_CAL_GNSS = {'BDS': 2.4, 'GAL': 2.4 ,'GLO': 3.8, 'GPS': 2.7}

# Uncertainty arising from choice of UTC model
# Values from Defraigne et al 2023 Metrologia 60 065010,  DOI : 10.1088/1681-7575/ad0562
U_NAVMSG_GNSS = {'BDS': 0.2, 'GAL': 0.1 ,'GLO':1.2, 'GPS': 1.3}

CLOCK_5071_STD = 1 # standard model
CLOCK_5071_HPT = 2 # high performance model 
CLOCK_MASER    = 3

# Indices for computed data array
D_UTCK_BUTC = 0
U_UTCK_BUTC = 1
D_UTC_BUTC  = 2
U_UTC_BUTC  = 3

CIRT_REF_ISSUE = 421
CIRT_REF_MJD   = 59944

CIRT_WEB_API  = 0
CIRT_REPORT   = 1
CIRT_DOWNLOAD = 2

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

# Note that this does not assume that the input MJD is an integer
def MJDtoBDSWeekDay(mjd):
	# Epoch for BDT is  0h UTC 1 Jan 2006 == MJD 53736
	ttBDS = (mjd - 53736)*86400 # number of seconds since the epoch
	BDSWn = int(ttBDS/(7*86400))
	BDSday = int((ttBDS - 7*86400*BDSWn)/86400)
	return (BDSWn,BDSday)
	
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

# Note that leapSecs is added to to the time we calculate for
# ------------------------------------------
def DeltaGNSSUTC(gnss,Wn,Dn,leapSecs,tsCorr):
	deltaUTC = 0
	if (gnss == 'GPS') or (gnss == 'GAL') or (gnss == 'BDS'):
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

# Returns Circular T in a dictionary
# A current difficulty is that the Web api only returns the total uncertainty
# For most labs though uA << uB so uB = uTotal
# ---------------------------------------------
def LoadUTCk(lab,startMJD,stopMJD):
	# Code for testing - be kind to the BIPM web server
	ottp.Debug('Reading local file for UTCk data')
	data = {}
	fin = open(os.path.join(root,'reports/cirt.txt'),'r')
	firstMJD = lastMJD = None
	for l in fin:
		if l[0] == '#':
			continue
		vals = l.strip().split()
		if (len(vals) == 3):
			data[int(vals[0])] = [int(vals[0]),float(vals[1]),float(vals[2]),
				0,math.sqrt(float(vals[2])**2 - 0**2 )] #  TODO when the web API reports (uA,uB)
			if not firstMJD:
				firstMJD = int(vals[0])
			lastMJD = int(vals[0])
	fin.close()
	return data,firstMJD,lastMJD

# ---------------------------------------------
def FetchUTCk(lab,startMJD,stopMJD):
	ottp.Debug(f'Fetching Circular T data for {lab}');
	ottp.Debug(f'{httpRequest}scale=utc&lab={lab}&mjd1={startMJD}&mjd2={stopMJD}&outfile=txt')
	try:
		r = requests.get(f'{httpRequest}scale=utc&lab={lab}&mjd1={startMJD}&mjd2={stopMJD}&outfile=txt',verify=rootCert)
	except:
		return None,None,None
	
	lines = r.text.split('\r\n')
	data = {}
	firstMJD = lastMJD = None
	for l in lines:
		l = l.strip()
		if not(l):
			continue
		if l[0] == '#':
			continue
		vals = l.strip().split() 
		# If there is no data for an MJD, then nothing is returned
		# If the uncertainty is unavailable, 0 is returned
		if (len(vals) == 3): # Three values are always returned
			if vals[2] == '0.0':
				data[int(vals[0])] = [int(vals[0]),float(vals[1]),None,None,None]
			else:	
				data[int(vals[0])] = [int(vals[0]),float(vals[1]),float(vals[2]),
					0,math.sqrt(float(vals[2])**2 - 0**2 )] #  TODO when the web API reports (uA,uB)
			if not firstMJD:
				firstMJD = int(vals[0])
			lastMJD = int(vals[0])
	return data,firstMJD,lastMJD

# ---------------------------------------------
def LoadCircularT(lab,startMJD,stopMJD):
	ottp.Debug('Reading Circular T from local mirror')
	# Estimate the range of issues to load
	startIssue = CIRT_REF_ISSUE + int(math.floor((startMJD - CIRT_REF_MJD)/(365.25/12.0)))
	stopIssue  = CIRT_REF_ISSUE + int(math.ceil((stopMJD - CIRT_REF_MJD)/(365.25/12.0))) # may go one too far but that's OK
	ottp.Debug(f'{startMJD},{stopMJD}->issues {startIssue},{stopIssue}')
	utck = {} # using a dictionary gets rid of the duplicates automatically
	labregex = r'^' + lab
	firstMJD = lastMJD = None
	for issue in range(startIssue,stopIssue+1):
		fName = os.path.join(cirtDir,f'cirt.{issue}')
		if os.path.exists(fName):
			ottp.Debug(f'Reading {fName}')
			fin = open(fName,'r')
			mjds = []
			for l in fin:
				if (re.search(r'uA\s+uB\s+u\s*$',l)):
					args = l.strip().split()
					for i in range(1,len(args)-3):
						#utck[int(args[i])] = [int(args[i]),None,None,None,None]
						mjds.append(int(args[i]))
					break
			if firstMJD == None:
				firstMJD = mjds[0]
			lastMJD = mjds[-1]
			for l in fin:
				if (re.search(labregex,l)):
					# Remove lab name and location for ease of parsing
					indx = l.find(')')
					ll = l[indx+1:]
					# Occasionally there are notes at the end, chop these off 
					# to make parsing easier
					indx = ll.find('(')
					if indx:
						ll = ll[0:indx]
					cols = ll.strip().split()
					utckm = []
					for i in range(0,len(mjds)):
						if (cols[i] == '-'): # missing data
							utckm.append(None) # tag for cleanup
						else:
							utckm.append(float(cols[i]))
					
					# If there is no data then no uncertainties are given in Circular T, logically enough
					if len(cols) == len(mjds): # no uncertainties are present
						uncerts = [None,None,None]
					else:
						uncerts = []
						for i in range(len(cols)-3,len(cols)):
							if (cols[i] == '-' or cols[i] == 'NC'): # missing/uncalibrated data
								uncerts.append(None)
							else:
								uncerts.append(float(cols[i]))
					for mi in range(len(mjds)):
						utck[mjds[mi]] = [ mjds[mi], utckm[mi], uncerts[2], uncerts[0],uncerts[1] ]  # note ordering of uncertainties
					break
			fin.close()
	
	return utck,firstMJD,lastMJD

		
# ---------------------------------------------
# This uses the list representation!
# 0->MJD, 1->UTC-UTCk, 2->u, 3->uA, 4->uB
def GetNearestCirtU(cirt,mjd):
	for i in range(0,len(cirt)-1):
		if (cirt[i][0] >= mjd and mjd <= cirt[i+1][0]):
			ottp.Debug(f'{cirt[i][0]} {cirt[i+1][0]} {mjd} {cirt[i][2]} ')
			return cirt[i][3], cirt[i][4]
	ottp.Debug(f'GetNearestCirtU failed for {mjd}: last MJD = {cirt[-1][0]}, using {cirt[-1][3]}  {cirt[-1][4]}')
	return cirt[-1][3],cirt[-1][4] # not available, use the last known values

# ---------------------------------------------
# Linear interpolation  of Circular T
# Returns the interpolated value and its uncertainty
# This uses the dict representation!

def InterpolateUTC(mjd,mjd0,mjd1,cirt):
	utc0 = cirt[mjd0][1]
	utc1 = cirt[mjd1][1]
	
	# Check the boundaries	
	if utc0 == None:
		if utc1 == None:
			ottp.Debug(f'InterpolateUTC: no data for {utc0} and {utc1}') # UNTESTED 
			return None,None 
		elif mjd == mjd1:
			ottp.Debug(f'InterpolateUTC: no data for {utc0} but MJD = {mjd1}, return {cirt[mjd1]} +/- 0')
			return cirt[mjd1][1],0 # no interpolation uncertainty
		else:
			ottp.Debug(f'InterpolateUTC: interpolation failure')
			return None,None
			
	if utc1 == None:
		if mjd == mjd0:
			ottp.Debug(f'InterpolateUTC: no data for {utc1} but MJD = {mjd0}, return {cirt[mjd1]} +/- 0') # UNTESTED 
			return cirt[mjd0][1],0 # no interpolation uncertainty
		else:
			ottp.Debug(f'InterpolateUTC: interpolation failure')
			return None,None
			
	dmjd = mjd - mjd0
	
	# The uncertainty is calculated by
	# using the recommended uncertainty for the mean frequency, as given in CCTF WGMRA Guideline 4 
	uFreq = cirt[mjd0][3] * SQRT2/ 5.0  # fractional, ns/day
	minD = dmjd # use the minimum distance 
	if mjd1-mjd < minD:
		minD = mjd1-mjd
	utcm = utc0 + dmjd*(utc1-utc0)/5 # estimate of UTC - UTC(k) by linear interpolation
	uc   = math.sqrt((uFreq*minD)**2 ) 
	ottp.Debug(f'InterpolateUTC {mjd} {utcm} +/- {uc}, [{utc0} {utc1}]')
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

configFile = os.path.join(root,'etc/butc.conf')
rnxDir = os.path.join(root,'rinex')
tmpDir  = os.path.join(root,'tmp')
cggttsDir = os.path.join(root,'cggtts')
cirtDir = os.path.join(root,'cirt')
db = os.path.join(root,'butcgnss.db')
cirtSource = CIRT_REPORT 
nPrevDays = NPREVDAYS

lab = 'AUS'
httpRequest = 'https://webtai.bipm.org/api/v1.0/get-data.html?'
rootCert = None # if you use an empty string, this will skip SSL verification which is a bad thing

winSize = REFSYS_AVG_WINDOW

clockModel = CLOCK_5071_STD

if ottp.LibMajorVersion() >= 0 and ottp.LibMinorVersion() < 2: 
	sys.exit('Need ottplib minor version >= 2')

examples='TO DO'
parser = argparse.ArgumentParser(description='Generate UTC(k) - bUTC_GNSS and UTC - bUTC_GNSS ',
	formatter_class=argparse.RawDescriptionHelpFormatter,epilog=examples)
parser.add_argument('mjd',nargs = '*',help='first MJD [last MJD] (if not given, the MJD of the previous day is used as the last MJD)')
parser.add_argument('--config','-c',help='use an alternate configuration file',default=configFile)
parser.add_argument('--debug','-d',help='debug (to stderr)',action='store_true')
parser.add_argument('--cirt',help='source of data for Circular T (webapi,report,download)')
parser.add_argument('--nprev',help='number of previous days to process',default=nPrevDays)
parser.add_argument('--version','-v',help='show version and exit',action='store_true')

args = parser.parse_args()

if (args.version):
	ShowVersion()
	exit()

if (args.config):
	configFile = args.config
	if (not os.path.isfile(configFile)):
		ottp.ErrorExit(configFile + ' not found')

if (args.cirt):
	if args.cirt == 'webapi':
		cirtSource = CIRT_WEB_API
	elif args.cirt == 'report':
		cirtSource = CIRT_REPORT
	elif args.cirt == 'download':
		cirtSource = CIRT_DOWNLOAD
	else:
		ottp.ErrorExit(f'Bad option --cirt {args.cirt}')

nPrevDays = int(args.nprev)

debug = args.debug
ottp.SetDebugging(debug)
cggtts.SetWarnings(debug)

# If an MJD range is not specified then
# the processed range is the previous day
# and 60 or so days before that to pick up Circular T

mjdToday  = ottp.MJD(time.time())
stopMJD   = mjdToday - 1 # previous day
startMJD  = stopMJD - nPrevDays

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

cfg=ottp.Initialise(configFile,['main:gnss'])

gnss = cfg['main:gnss'].split(',')
gnss = [g.strip() for g in gnss] 

if 'main:root' in cfg:
	root = cfg['main:root']

if 'main:lab' in cfg:
	lab = cfg['main:lab']
	
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
	
if ('main:root certificate' in cfg):
	rootCert= cfg['main:root certificate']

if ('database:file' in cfg):
	db = ottp.MakeAbsoluteFilePath(cfg['database:file'],root,os.path.join(home,'database'))

if 'rinex:path' in cfg:
	rnxDir = ottp.MakeAbsolutePath(cfg['rinex:path'],root)
	# check other necessary things have been defined
	if not(CheckConfig(cfg,['rinex:station name','rinex:version'])):
		ottp.ErrorExit('Missing entries in configuration file')
	staName = cfg['rinex:station name']
	rnxVersion = int(cfg['rinex:version'])

# Connect to the database
# Do this before we download CircularT because we need uA from the database

ottp.Debug(f'Connecting to the database {db}')
dbc = sqlite3.connect(db) # this opens a connection to the database, creating the db if it doesn't exist
curs = dbc.cursor()       #

# Create the table, if it doesn't exist
# The primary key for the table is the MJD
# This has four entries for each GNSS, all kept as unrounded REALs
# Rounding etc is done when reports are created 
# UTCk - bUTC_GNSS , u, UTC - buTC_gnss,  u
# Also keep track of the values of uA and uB used
# noting that the WebAPI currently only returns the quadrature sum

curs.execute(
		"CREATE TABLE IF NOT EXISTS butcgnss ("
    "MJD INTEGER PRIMARY KEY,"
    "UTCk_BDS REAL,UTCk_BDS_u REAL,UTC_BDS REAL,UTC_BDS_u REAL,"
    "UTCk_GAL REAL,UTCk_GAL_u REAL,UTC_GAL REAL,UTC_GAL_u REAL,"
    "UTCk_GLO REAL,UTCk_GLO_u REAL,UTC_GLO REAL,UTC_GLO_u REAL,"
    "UTCk_GPS REAL,UTCk_GPS_u REAL,UTC_GPS REAL,UTC_GPS_u REAL)"
)


ottp.Debug(f'Processing range is {startMJD} - {stopMJD}')

# CircularT data is downloaded each run
# One year of data is a 2K file so we shouldn't worry too much about 
# about the load on the web server
# Need to get Circular T so that we know the current uncertainty uB
# Unfortunately the Web API returns the total uncertainty.
# In most cases, and for us, uB dominates so for the moment
# we'll use the total uncertainty for uB. Typically, we won't have a matching uB
# for an MJD, so we'll just be using the most recent anyway (noting that this will be fixed
# later when the UTC update is done)

# Initial guess on the range of Circular T to ask for is
# [startMJD -7, stopMJD + 7] say

cirtStartMJD = startMJD - 7
cirtStopMJD  = stopMJD  + 7

# FIXME when an MJD range is specified, this may need to be extended backwards 
# to get a reported day

ottp.Debug(f'Getting Circular T data for {cirtStartMJD} - {cirtStopMJD}')

if cirtSource == CIRT_WEB_API:
	cirt, firstMJD,lastMJD = FetchUTCk(lab,cirtStartMJD,cirtStopMJD)	
elif cirtSource == CIRT_REPORT:
	cirt, firstMJD,lastMJD = LoadCircularT(lab,cirtStartMJD,cirtStopMJD)	
elif cirtSource == CIRT_DOWNLOAD:
	cirt,firstMJD,lastMJD = LoadUTCk(lab,cirtStartMJD,cirtStopMJD)

if not cirt:
	sys.exit("Couldn't get Circular T data")

cirtAsList = list(cirt.values())

prevCGGTTSFile = {}
prevCGGTTSFile['BDS'] = CGGTTS(None,None)
prevCGGTTSFile['GAL'] = CGGTTS(None,None)
prevCGGTTSFile['GLO'] = CGGTTS(None,None)
prevCGGTTSFile['GPS'] = CGGTTS(None,None)
newData = {}

mjd = startMJD - 2  # start two days earlier because of 1. immediate increment and 2. need the previous day's data

while mjd < stopMJD :
	mjd += 1 # doing this here means not having to increment multiple times
	
	ts = (mjd - 40587)*86400 # convert MJD to UNIX time
	dt = datetime.fromtimestamp(ts, tz=timezone.utc) # get date in UTC
	yyyy = dt.year
	mm   = dt.month
	dd   = dt.day
	doy = int(dt.strftime('%j'))
	
	ottp.Debug(f'\nProcessing {mjd}: {yyyy}-{mm}-{dd}')
	
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
	
	newData[mjd] = {} # we get an empty entry first time through but that's OK
	
	# 
	uAcirt,uBcirt = GetNearestCirtU(cirtAsList,mjd) # could be None
	
	for g in gnss:
		ottp.Debug(f'Processing {g}')
		
		if g == 'BDS':
			Wn,Dn = MJDtoBDSWeekDay(mjd)
		else:
			Wn,Dn = ottp.MJDtoGPSWeekDay(mjd)
		
		ottp.Debug(f'Processing {g} WN = {Wn} DN = {Dn}')
	
		newData[mjd][g] = [None,None,None,None] # Fields are UTC(k)-bUTC_GNSS, u, UTC - bUTC_GNSS, u]
		
		fName = FindCGTTSFile(g,mjd)
		
		if fName:
			ottp.Debug(f'Reading {fName}')
			cgf = CGGTTS(fName,mjd)
			cgf.Read()
		else:
			ottp.Debug('No CGGTTS file found')
			cgf = CGGTTS(None,mjd)
		
		# The first one is for the day before the start MJD so no more to do
		if (mjd==startMJD-1):
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
			
		deltaUTC = DeltaGNSSUTC(g,Wn,Dn,leapSecs,tsCorr)
			
		newData[mjd][g][D_UTCK_BUTC] = refsys0 + deltaUTC
		if uBcirt == None:
			newData[mjd][g][U_UTCK_BUTC] = None # UNTESTED
			ottp.Debug('UTCk - bUTC {} {:g} {:g} +/- ? refsys0={:g} urefsys0={:g}'.format(g,mjd, newData[mjd][g][D_UTCK_BUTC],refsys0,uRefsys0))
		else:
			newData[mjd][g][U_UTCK_BUTC] = math.sqrt(uRefsys0**2 + uBcirt**2 + U_CAL_GNSS[g]**2 +  U_NAVMSG_GNSS[g]**2) # only uB is relevant here
			ottp.Debug('UTCk - bUTC {} {:g} {:g} +/- {:g} refsys0={:g} urefsys0={:g}'.format(g,mjd, newData[mjd][g][D_UTCK_BUTC],newData[mjd][g][U_UTCK_BUTC],refsys0,uRefsys0))
		
		# Now update UTC - bUTC_GNSS, if we can
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
				newData[mjd][g][D_UTC_BUTC] = cirt[mjd][1] + newData[mjd][g][D_UTCK_BUTC]
				if uAcirt == None or uBcirt == None:
					newData[mjd][g][U_UTC_BUTC] = None # untested
				else:
					newData[mjd][g][U_UTC_BUTC] = math.sqrt(uRefsys0**2 +  uAcirt**2 + uBcirt**2 + U_CAL_GNSS[g]**2 +  U_NAVMSG_GNSS[g]**2) # no contribution from instability
				mjdLastupdateUTC = mjd
			prevCGGTTSFile[g] = cgf # this is correct - we only have mjd0, so we're done after tewsing for mjd==mjd0 
			continue
		
		utc0 = cirt[mjd0][1]
		utc1 = cirt[mjd1][1]	
		
		utcDiff,utcInterpUncert = InterpolateUTC(mjd,mjd0,mjd1,cirt)
		if utcDiff == None:
			ottp.Debug(f'UTC interpolation failed on [{mjd0},{mjd1}]')
			newData[mjd][g][D_UTC_BUTC] = None
		else:
			newData[mjd][g][D_UTC_BUTC] = utcDiff + newData[mjd][g][D_UTCK_BUTC]
			# Uncertainty sources:
			# time transfer noise == uRefsys0
			# UTC interpolation uncertainty == utcUncert (which includes the link calibration uncertainty)
			# UTC(k) instability 
			# GNSS provider link's calibration uncertainty
			# UTC prediction 
			UTCkInstability = ClockStability(clockModel,mjd - mjd0)
			if  ClockStability(clockModel,mjd - mjd0) <  UTCkInstability:
				UTCkInstability = ClockStability(clockModel,mjd1- mjd)
			
			newData[mjd][g][U_UTC_BUTC] = math.sqrt(uRefsys0**2 + cirt[mjd0][3]**2 + cirt[mjd0][4]**2 + utcInterpUncert**2 +
					UTCkInstability**2+ U_CAL_GNSS[g]**2 +  U_NAVMSG_GNSS[g]**2) 
			
		prevCGGTTSFile[g] = cgf # save it for the next iteration
	
	if (mjd == startMJD):
		continue

#print(f"SQLite Library Version: {sqlite3.sqlite_version}")
#print(f"Python sqlite3 Module Version: {sqlite3.version}")

# Finally, update the database
for m in newData:
	# Build an 'upsert' command
	inscmd = 'INSERT INTO butcgnss (mjd'
	valscmd = f'VALUES ({m}'
	updatedcols = ''
	cnt = 0
	for g in gnss:
		inscmd += f',UTCk_{g},UTCk_{g}_u,UTC_{g},UTC_{g}_u'
		for i in range(0,4):
			if newData[m][g][i] == None:
				valscmd += ',NULL'
			else:
				valscmd += f',{newData[m][g][i]}'
		if cnt>0:
			updatedcols += f',\nUTCk_{g}=EXCLUDED.UTCk_{g},'
		else:
			updatedcols += f'\nUTCk_{g}=EXCLUDED.UTCk_{g},'
		updatedcols += f'\nUTCk_{g}_u=EXCLUDED.UTCk_{g}_u,'
		updatedcols += f'\nUTC_{g}=EXCLUDED.UTC_{g},'
		updatedcols += f'\nUTC_{g}_u=EXCLUDED.UTC_{g}_u'
		cnt += 1
	
	inscmd += ')\n'
	valscmd += ')\n'

	cmd = inscmd + valscmd + 'ON CONFLICT(mjd)\n' + 'DO UPDATE SET' + updatedcols + ';\n' 
	curs.execute(cmd)
	
dbc.commit()
curs.close()
dbc.close()
  
ottp.Debug('Done!')
