#!/usr/bin/python

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
import shutil
import sys
import time

# This is where cggttslib is installed
sys.path.append("/usr/local/lib/python3.8/site-packages")  # Ubuntu 20.04
sys.path.append("/usr/local/lib/python3.10/site-packages") # Ubuntu 22.04

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
	
VERSION = '0.0.3'
AUTHORS = 'Michael Wouters'

REFSYS_AVG_WINDOW = 2 # window size for averaging

U_GNSS_LINK = 2.5     # uncertainty of GNSS provider's link calibration 
U_BUTC_MODEL = 1      # uncertainty arising from choice of UTC model

# ------------------------------------------
def ShowVersion():
	print (os.path.basename(sys.argv[0])+" "+VERSION)
	print ('Written by ' + AUTHORS)
	return

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
# ------------------------------------------
def GetTimeSysCorr(gnss,navFile):
	
	fnav = open(navFile,'r')
	for l in fnav:
		
		if gnss == 'GPS':
			if l[0:4] == 'GPUT': # covers version 2.12 ->
				# Fields are offset a0, rate a1, tow, wn
				return [float(l[5:22]), float(l[22:38]), int(l[38:45]),int(l[45:50])]
			# TODO earlier RINEX 
				
	fnav.close()

# -------------------------------------------
# TDEV for Cs  5071
# 
def TDEV_5071_Std(tauDays):
	return 2.0*math.sqrt(tauDays) # in ns 
	
#----------------------------------------------
def CheckConfig(cfg,req):
	ok = True
	for r in req:
		if not(r in cfg):
			print(f'{r} is not set')
			ok = False
	return ok

# ---------------------------------------------
def GetRefsys(cgBef,cgAft,winSize):

	wSz = 1800*winSize 
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
	
	return refsys0,math.sqrt(urefsys0/(nTracks-1)),nTracks

# ---------------------------------------------
def formatNumber(fval,res=0.1):
	invres = 1/res
	return str(round(fval*invres)/invres) # note that this is banker's rounding
	
# ---------------------------------------------

#import requests
#r = requests.get("https://webtai.bipm.org/api/v1.0/get-data.html?scale=utc&lab=AOS")

# ---------------------------------------------
def GetCircularT(lab,startMJD,stopMJD):
	# Temporary code
	data = {}
	fin = open(os.path.join(root,'report/cirt.txt'),'r')
	firstMJD = lastMJD = None
	for l in fin:
		if l[0] == '#':
			continue
		vals = l.strip().split()
		if (len(vals) == 3):
			data[int(vals[0])] = [float(vals[1]),float(vals[2])]
			if not firstMJD:
				firstMJD = int(vals[0])
			lastMJD = int(vals[0])
	fin.close()
	return data,firstMJD,lastMJD
	
# ---------------------------------------------
def WriteHeader(fout):
	fout.write('MJD\n')

# ---------------------------------------------
def WriteFooter(fout,footer):
	if (os.path.isfile(footer)):
		with open(footer, 'r') as fin:
			txt = fin.read()
			fout.write('\n\n')
			fout.write(txt)
	else:
		ottp.Debug('Unable to open {footer}')

# ---------------------------------------------
# Linear interpolation  of Circular T
# Returns the interpolated value and its uncertainty
def InterpolateUTC(mjd,mjd0,mjd1,cirt):
	utc0 = cirt[mjd0][0]
	utc1 = cirt[mjd1][0]	
	dmjd = mjd - mjd0
	return utc0 + dmjd*(utc1-utc0)/5,math.sqrt((cirt[mjd0][1]*(1-dmjd/5))**2 + (cirt[mjd1][1]*dmjd/5)**2)
	
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

winSize = REFSYS_AVG_WINDOW

minUTCkUncertainty = 5
minUTCUncertainty = 6

if ottp.LibMajorVersion() >= 0 and ottp.LibMinorVersion() < 2: 
	sys.exit('Need ottplib minor version >= 2')

examples='TO DO'
parser = argparse.ArgumentParser(description='Generate UTC(k) - bUTC_GNSS and UTC - bUTC_GNSS ',
	formatter_class=argparse.RawDescriptionHelpFormatter,epilog=examples)
parser.add_argument('mjd',nargs = '*',help='first MJD [last MJD] (if not given, the MJD of the previous day is used)')
parser.add_argument('--config','-c',help='use an alternate configuration file',default=configFile)
parser.add_argument('--utc','-u',help='generate UTC-bUTC_gnss differences for the previous month',action='store_true')
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
cggtts.SetWarnings(debug)

cfg=ottp.Initialise(configFile,['main:gnss']);

gnss = cfg['main:gnss'].split(',')

if 'main:root' in cfg:
	root = cfg['main:root']

if 'main:report path' in cfg:
	reportDir = ottp.MakeAbsolutePath(cfg['main:report path'],root)
	
if 'main:footer' in cfg:
	footer = ottp.MakeAbsoluteFilePath(cfg['main:footer'],root,footer)

if 'main:utck uncertainty' in cfg:
	minUTCkUncertainty = float(cfg['main:utck uncertainty'])
	
if 'main:utc uncertainty' in cfg:
	minUTCUncertainty = float(cfg['main:utc uncertainty'])
	
if 'rinex:path' in cfg:
	rnxDir = ottp.MakeAbsolutePath(cfg['rinex:path'],root)
	# check other necessary things have been defined
	if not(CheckConfig(cfg,['rinex:station name','rinex:version'])):
		ottp.ErrorExit('Missing entries in configuration file')
	staName = cfg['rinex:station name']
	rnxVersion = int(cfg['rinex:version'])

dailyUpdate = False
UTCupdate   = False
if args.utc:
	if args.mjd:
		sys.exit("An MJD range cannot be used with the UTC option")
	
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
	cirt,firstMJD,lastMJD = GetCircularT('AUS',startMJD - 7,stopMJD + 7) 
	if not cirt:
		sys.exit("Couldn't get Circular T")

	# Check that the required MJDs are there
	if not(firstMJD <= startMJD and stopMJD - 5  <= lastMJD ): # probably won't get the last day of the month
		sys.exit(f"Not all required MJDs are available in Circular T (wanted [{startMJD},{stopMJD} - 5], got [{firstMJD},{lastMJD}])")
	
	UTCupdate = True
	
else:
	ts = time.time() - 86400 # previous day
	# Set the default processing range
	stopMJD = ottp.MJD(ts)  
	
	# Now need to get
	dt = datetime.fromtimestamp(ts, tz=timezone.utc) # get date in UTC
	ts   = ts - 86400*(dt.day-1) # FIXME check month rollover
	startMJD = ottp.MJD(ts)-1 # need to get the previous day as well
	
	dailyUpdate = True 
	ottp.Debug(f'Default processing range {startMJD} - {stopMJD}')

	if (args.mjd):
		if 1 == len(args.mjd):
			startMJD = int(args.mjd[0]) - 1 # need to get the previous day 
			stopMJD  = int(args.mjd[0])        # 
		elif ( 2 == len(args.mjd)):
			startMJD = int(args.mjd[0]) - 1 # need to get the previous day 
			stopMJD  = int(args.mjd[1])
			if (stopMJD < startMJD):
				ottp.ErrorExit('Stop MJD is before start MJD')
		else:
			ottp.ErrorExit('Too many MJDs')
		dailyUpdate = False


currReport = None # initially, the report file is closed or doesn't yet exist
mjdLastUTCupdate = None # for tracking if we updated UTC - buTC_GNSS for the whople month
prevCGGTTSFile = {}
prevCGGTTSFile['GPS'] = CGGTTS(None,None)
newData = {}

for mjd in range(startMJD,stopMJD + 1):
	
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
				currReport.write('UTC - bUTC_GNSS updated\n')
			currReport.write('####################################################################\n')
			currReport.close()
			
		reportName = os.path.join(reportDir,f'brutc{mm:02d}{yyyy:4d}.txt')
		ottp.Debug(f'Creating report {reportName}')
		try:
			currReport = open(reportName,"w")
		except:
			sys.exit(f"Couldn't open {reportName}")
		if UTCupdate:
			# Check whether the file needs processing
			pass
		WriteHeader(currReport)
		
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
			refsys0,uRefsys0,nTracks = GetRefsys(prevCGGTTSFile[g],cgf,winSize);
			#print(refsys0,uRefsys0,nTracks)
		
		if (refsys0 == None):
			ottp.Debug('Insufficent CGGTTS data')
			prevCGGTTSFile[g] = cgf
			continue
	
		# Find the navigation file
		navFile,zExt = rinex.FindNavigationFile(rnxDir,staName,yyyy,doy,rnxVersion,False) # don't exit if not found

		if not navFile:
			ottp.Debug("Couldn't find a navigation file for {yyyy} {doy}")
			prevCGGTTSFile[g] = cgf
			continue
			
		ottp.Debug(f'Found {navFile}, compression = {zExt}')
		src = navFile + zExt
		srcBase = os.path.basename(src)
		dst = os.path.join(tmpDir,srcBase)
		# We may not own the file so we need to make a local copy
		shutil.copy(src,dst)
		navFile,zAlgorithm = rinex.Decompress(dst)
		tsCorr = GetTimeSysCorr(g,navFile)
		leapSecs =rinex.GetLeapSeconds(navFile,rnxVersion)
		if g == 'GPS':
			# From the ICD
			# t_UTC = t_E - delta_UTC where t_E is 'effective' GPS time 
			# delta_UTC = dt_LS + A0 + A1*(t_E - t0t + 604800*(WN - Wn_t) )
			# We want the value at UTC0 for the day which means
			# t_E = 86400*gpsDn + leapSecs
			deltaUTC = tsCorr[0] + tsCorr[1]*(86400*gpsDn + leapSecs - tsCorr[2] + 604800*(gpsWn - tsCorr[3]))
			#print(deltaUTC,tsCorr,gpsDn,gpsWn,leapSecs)
			reportData[g][0] = refsys0 + deltaUTC*1.0E9
			reportData[g][1] = math.sqrt(uRefsys0**2 + 9 + U_GNSS_LINK**2 +  U_BUTC_MODEL**2) # TODO link uncertainty
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
			if not(mjd1 in cirt):
				ottp.Debug(f'{mjd1} not in CirT')
				if (mjd == mjd0):
					reportData[g][2] = cirt[mjd][0] + reportData[g][0]
					reportData[g][3] = math.sqrt(uRefsys0**2 + cirt[mjd][1]**2 + U_GNSS_LINK**2 +  U_BUTC_MODEL**2) # no contribution from instability
					mjdLastUTCupdate = mjd
				prevCGGTTSFile[g] = cgf
				continue
				
			utc0 = cirt[mjd0][0]
			utc1 = cirt[mjd1][0]	
			
			utcDiff,utcUncert = InterpolateUTC(mjd,mjd0,mjd1,cirt)
			reportData[g][2] = utcDiff + reportData[g][0]
			# Uncertainty sources:
			# time transfer noise == uRefsys0
			# UTC interpolation uncertainty == utcUncert (which includes the link calibration uncertainty)
			# UTC(k) instability 
			# GNSS provider link's calibration uncertainty
			# UTC prediction 
			UTCkInstability = TDEV_5071_Std(mjd - mjd0); 
			if  TDEV_5071_Std(mjd1- mjd) <  UTCkInstability:
				UTCkInstability = TDEV_5071_Std(mjd1- mjd)
			
			#print(mjd,uRefsys0, utcUncert,UTCkInstability,U_GNSS_LINK,U_BUTC_MODEL)
			reportData[g][3] = math.sqrt(uRefsys0**2 + utcUncert**2 + UTCkInstability**2+ U_GNSS_LINK**2 +  U_BUTC_MODEL**2) # FIXME
			
			mjdLastUTCupdate = mjd
		prevCGGTTSFile[g] = cgf # save it for the next time
	
	if (mjd == startMJD):
		continue
	
	# If we can't generate the data, mark as missing'
	outputLine = f'{yyyy:4d}-{mm:02d}-{dd:02d} {mjd:5d}'
	missingData = '*'
	for g in gnss:
		if reportData[g][0]==None:
			outputLine += f'{missingData:>9}{missingData:>5}{missingData:>9}{missingData:>5}'
		elif reportData[g][2]==None:
			reportedUTCkUncert = reportData[g][1]
			if reportedUTCkUncert  < minUTCkUncertainty:
				reportedUTCkUncert = minUTCkUncertainty
			outputLine += '{:>9}{:>5}{:>9}{:>5}'.format(formatNumber(reportData[g][0]),formatNumber(reportedUTCkUncert,1),missingData,missingData)
		else:
			reportedUTCkUncert = reportData[g][1]
			if reportedUTCkUncert  < minUTCkUncertainty:
				reportedUTCkUncert = minUTCkUncertainty
			reportedUTCUncert = reportData[g][3]
			if reportedUTCUncert  < minUTCUncertainty:
				reportedUTCUncert = minUTCUncertainty
			outputLine += '{:>9}{:>5}{:>9}{:>5}'.format(formatNumber(reportData[g][0],1),formatNumber(reportedUTCkUncert,1),formatNumber(reportData[g][2]),formatNumber(reportedUTCUncert,1))
	if (dailyUpdate or UTCupdate):
		currReport.write(outputLine+'\n')
	newData[mjd] = outputLine + '\n'

# Finish off any currently open file	
if currReport:
	# TODO want to add a note to tag that UTC - bUTC_GNSS is updated
	WriteFooter(currReport,footer)
	currReport.close()		

	
	
	
