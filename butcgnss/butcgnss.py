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
	
VERSION = '0.0.2'
AUTHORS = 'Michael Wouters'

REFSYS_AVG_WINDOW = 2 # window size for averaging

# ------------------------------------------
def ShowVersion():
	print (os.path.basename(sys.argv[0])+" "+VERSION)
	print ('Written by ' + AUTHORS)
	return

# ------------------------------------------
# Writes header of the report file
def WriteReportHeader(fout):
	pass

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
def formatNumber(fval):
	return str(round(fval*10)/10.0) # in 0.1 ns note that this is banker's rounding
	

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

winSize = REFSYS_AVG_WINDOW

if ottp.LibMajorVersion() >= 0 and ottp.LibMinorVersion() < 2: 
	print('ottplib minor version < 2')
	sys.exit(1)

examples='TO DO'
parser = argparse.ArgumentParser(description='Generate UTC(k) - bUTC_GNSS and UTC - bUTC_GNSS ',
	formatter_class=argparse.RawDescriptionHelpFormatter,epilog=examples)
parser.add_argument('mjd',nargs = '*',help='first MJD [last MJD] (if not given, the MJD of the previous day is used)')
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
cggtts.SetWarnings(debug)

cfg=ottp.Initialise(configFile,['main:gnss']);

gnss = cfg['main:gnss'].split(',')

if 'main:root' in cfg:
	root = cfg['main:root']

if 'main:report path' in cfg:
	reportDir = ottp.MakeAbsolutePath(cfg['main:report path'],root)
	
if 'rinex:path' in cfg:
	rnxDir = ottp.MakeAbsolutePath(cfg['rinex:path'],root)
	# check other necessary things have been defined
	if not(CheckConfig(cfg,['rinex:station name','rinex:version'])):
		ottp.ErrorExit('Missing entries in configuration file')
	staName = cfg['rinex:station name']
	rnxVersion = int(cfg['rinex:version'])
	
ts = time.time() - 86400 # previous day
# Set the default processing range
stopMJD = ottp.MJD(ts) 
dt = datetime.fromtimestamp(ts, tz=timezone.utc) # get date in UTC
ts   = ts - 86400*(dt.day-1) # FIXME check month rollover
startMJD = ottp.MJD(ts)-1 # need to get the previous day 
auto = True # i
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
	auto = False # if we manually specify a range then the existing data are overwritten


# First pass - update UTC(k) - bUTC_GNSS
freport = None # initially, the report file is closed or doesn't yet exist
prevCGGTTSFile = {}
prevCGGTTSFile['GPS'] = CGGTTS(None,None)

for mjd in range(startMJD,stopMJD + 1):
	
	reportData = {}
	
	ts = (mjd - 40587)*86400 # convert MJD to UNIX time
	dt = datetime.fromtimestamp(ts, tz=timezone.utc) # get date in UTC
	yyyy = dt.year
	mm   = dt.month
	dd   = dt.day
	doy = int(dt.strftime('%j'))
	gpsWn,gpsDn = ottp.MJDtoGPSWeekDay(mjd)
	ottp.Debug(f'\nProcessing {mjd}: {yyyy}-{mm}-{dd}, WN={gpsWn} DN={gpsDn}')
	
	# Is the required monthly report open?
	if freport == None:
		# First, does it exist ?
		reportName = os.path.join(reportDir,f'brutc{mm:02d}{yyyy:4d}.txt')
		if (not os.path.isfile(reportName)):
			ottp.Debug(f'Creating new report {reportName}')
			freport = open(reportName,'w')
		else:
			pass
	
	for g in gnss:
		ottp.Debug(f'Processing {g}')
		reportData[g] = [None,None,None,None] # Fields are UTC(k)-bUTC_GNSS, u, UTC - bUTC_GNSS, u]
		
		fname = FindCGTTSFile(g,mjd)
		
		if fname:
			ottp.Debug(f'Reading {fname}')
			cgf = CGGTTS(fname,mjd)
			cgf.Read()
		else:
			ottp.Debug('No CGGTTS file found')
			cgf = CGGTTS(None,mjd)
		
		# The first one is for the day before the start MJD so no more to do
		if (mjd==startMJD):
			prevCGGTTSFile[g] = cgf # save it 
			continue
		
		# Calculate the average REFSYS, if we can
		if prevCGGTTSFile[g].fileName or cgf.fileName:
			refsys0,uRefsys0,nTracks = GetRefsys(prevCGGTTSFile[g],cgf,winSize);
			print(refsys0,uRefsys0,nTracks)
		
			if not (refsys0 == None):
				# Find the navigation file
				navFile,zExt = rinex.FindNavigationFile(rnxDir,staName,yyyy,doy,rnxVersion,False) # don't exit if not found

				if navFile:
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
						# print(deltaUTC,tsCorr,gpsDn,gpsWn,leapSecs)
						reportData[g][0] = refsys0 + deltaUTC*1.0E9
						reportData[g][1] = uRefsys0
		else:
			ottp.Debug('Insufficent CGGTTS data')
	
		prevCGGTTSFile[g] = cgf # save it for the next time
		
	# If we can't generate the data, mark as missing'
	outputLine = f'{mjd:5d}'
	missingData = '*'
	for g in gnss:
		if reportData[g][0]==None:
			outputLine += f'{missingData:^9}{missingData:^9}'
		else:
			outputLine += '{:>9}{:>9}'.format(formatNumber(reportData[g][0]),formatNumber(reportData[g][1]))
	print(outputLine)
# Second pass - update UTC - bUTC_GNSS
