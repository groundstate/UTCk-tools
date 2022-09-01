#!/usr/bin/python3
#

#
# The MIT License (MIT)
#
# Copyright (c) 2022 Michael J. Wouters
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
# Report on planned steering of UTC(AUS) using Rapid UTC
# 



import allantools
import argparse
import datetime
import glob
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mticker
import os
import re
import requests
import socket
import sys
import time

import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.image import MIMEImage

# This is where ottplib is installed
sys.path.append("/usr/local/lib/python3.6/site-packages") # Ubuntu 18.04
sys.path.append("/usr/local/lib/python3.8/site-packages") # Ubuntu 20.04

import ottplib

VERSION = "0.0.3"
AUTHORS = "Michael Wouters"

UTCR_LATENCY = 3

# ------------------------------------------
def Debug(msg):
	if (debug):
		sys.stderr.write(msg+'\n')
	return

# ------------------------------------------
def Log(logFile,msg):
	Debug(msg)
	try:
		flog = open(logFile,'a')
		flog.write('{} {}\n'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.gmtime()),msg))
		flog.close()
		flog.close()
	except:
		Debug('Unable to log message')
	return

# ------------------------------------------
def Warn(msg):
	if (not args.nowarn):
		sys.stderr.write('WARNING! '+ msg+'\n')
	return

# ------------------------------------------
def ErrorExit(msg):
	sys.stderr.write(msg+'\n')
	sys.exit(0)
	
# ------------------------------------------
def Initialise(configFile):
	cfg=ottplib.LoadConfig(configFile,{'tolower':True})
	if (cfg == None):
		ErrorExit("Error loading " + configFile)
		
	# Check for required arguments
	reqd = []
	for k in reqd:
		if (not k in cfg):
			ErrorExit("The required configuration entry " + k + " is undefined")
		
	return cfg
	
# ------------------------------------------
def FudgeData(mjds,utcr):
	utcclk = 'cs2269'
	clkdata = '/mnt/logger3/clocks/' + utcclk
	delay = 99.86 # delay to be added to align HROG - Cs2269
	i = 0
	for m in mjds:
		fname = '{}/{:d}.{}'.format(clkdata,m,utcclk)
		if (os.path.exists(fname)):
			fin = open(fname,'r')
			for l in fin:
				if re.match('#',l):
					continue
				[tod,rdg] = l.split()
				Debug('{:d} {} {}'.format(m,tod,rdg))
				utcr[i] = utcr[i] - float(rdg) + delay
				Debug('{:d} {} {} -> {:g}'.format(m,tod,rdg,utcr[i]))
				break
			fin.close()
		else:
			pass
		i += 1

# ------------------------------------------

appName = os.path.basename(sys.argv[0])
debug = False
scheduleSteer = True
forceSteer = False
minFitPts = 5 # minimum acceptable number of points for a frequency offset estimate
maxClockFOffset = 10.0 # ns/day
maxPhaseSlew = 10.0 # ns/day
maxLatency = 5 # 
utcrRefMJD = 56473 # Last MJD for UTCr_1326
freqGain = 0.5 # gain/weight for component of frequency adjustment due to mean frequency offset 
phaseGain = 0.5 # gain/weight for component of frequency adjustment due to current phase offset

UTCID = 'AUS'
enableSteerFile = 'ENABLE_STEERING' 
lastUTCrFile  = 'LAST_UTCR_DOWNLOAD' # file containing MJD of last UTCr download 

home =os.environ['HOME'] + '/'
user =os.environ['USER'] # remember to define this in the user's crontab
configFile = os.path.join(home,'etc/utcsteersched.conf')
repDir  = os.path.join(home,'utcsteer/reports')
logDir = os.path.join(home,'utcsteer/logs')
controlDir = os.path.join(home,'utcsteer/control')
scheduledDir = os.path.join(controlDir,'scheduled_steer')
processedDir = os.path.join(controlDir,'processed_steers')
historyLength = 90 # in days
tmpDir = os.path.join(home,'utcsteer/tmp')
recipients = 'Michael.Wouters@measurement.gov.au'
email = True
bipmurl = 'https://webtai.bipm.org/api/v0.2-beta'

parser = argparse.ArgumentParser(description='Report on UTCr(k) from Rapid UTC data',
	formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--config','-c',help='use this configuration file',default=configFile)
parser.add_argument('--debug','-d',help='debug (to stderr)',action='store_true')
parser.add_argument('--force','-f',help='force steering',action='store_true')
parser.add_argument('--show','-s',help='show plots',action='store_true')
parser.add_argument('--version','-v',action='version',version = appName + ' ' + VERSION + '\n' + 'Written by ' + AUTHORS)

args = parser.parse_args()

debug = args.debug
configFile = args.config
forceSteer = args.force

cfg = Initialise(configFile)

if ('main:history' in cfg):
	historyLength = int(cfg['main:history'])

if ('main:utc id' in cfg):
	UTCID = cfg['main:utc id']

if ('main:email recipients' in cfg):
	recipients = cfg['main:email recipients']
	
if ('main:reports' in cfg):
	repDir = os.path.join(home,cfg['main:reports']) 
	
logFile = os.path.join(logDir,'utcsteer.log') # this log will be common to several scripts
Log(logFile,'running')

tt = time.time()
mjdToday = int(tt/86400)+40587
Debug('MJD today is {:d}'.format(mjdToday))
dt = datetime.datetime.today()
dow = dt.weekday()
yyyy = dt.year
mm  = dt.month
dd = dt.day


# UTCr is published on Wednesday 
# There's a 3 day lag (using rapid orbits and clocks for CGGTTS tweaking ?)
# Calculate the expected MJD in UTCr
nWeeks = int(math.floor((mjdToday - utcrRefMJD)/7))
lastMJD = utcrRefMJD + 7 * nWeeks
Debug('Last MJD in UTCr should be ' + str(lastMJD))

# Sanity checks
# Should be running at least 3 days after 'lastMJD'
if not(args.force):
	if (mjdToday - lastMJD < UTCR_LATENCY):
		ErrorExit('Running too early - next run needs to be at least {:d}'.format(lastMJD+3))
	# Shouldn't run too late either
	if (mjdToday - lastMJD > UTCR_LATENCY + 1):
		ErrorExit('Running too late (limit is MJD {:d}) '.format(lastMJD + UTCR_LATENCY + 1))
	
mjd1 = lastMJD - historyLength
mjd2 = lastMJD

GETDATA = True

# Get the data

dmjd = []
dutck = []

if GETDATA: # TEMPORARY
	Debug('Fetching data for interval {:d} {:d}'.format(mjd1,mjd2))
	
	# Have we already got the data ?
	# Check the file lastUTCrDownload
	fin = open(os.path.join(controlDir,lastUTCrFile),'r')
	lastUTCr = -1 # flags failure to get this
	for l in fin:
		if re.match(r'^#',l): # ignore comments
			continue
		matches = re.match(r'MJD\s+(\d{5})',l)
		if matches:
			lastUTCr = int(matches.group(1))
			Debug('Last UTCr {:d}'.format(lastUTCr))
			break
	fin.close()
	
	# We have already checked the run window (unless --force is in effect)
	if (lastUTCr == lastMJD):
		Debug('Nothing to do')
		sys.exit(0)
	
	httpreq = '{}/get-data.html?scale=utcr&lab={}&outfile=txt&mjd1={:d}&mjd2={:d}'.format(bipmurl,UTCID,mjd1,mjd2)
	try:
		resp = requests.get(httpreq)
	except:
		sys.exit(1)

	# Parse what we got back

	lines = resp.text.split('\r\n')
	for l in lines:
		if re.match('#',l): # ignore comments
			continue
		ldata = l.split()
		if len(ldata) == 2:
			try:
				dmjd.append(int(ldata[0]))
				dutck.append(float(ldata[1]))
			except:
				Debug('Bad data: ' + l)

	Debug('UTCr data: first = {:d}, last = {:d}'.format(dmjd[0],dmjd[-1]))

	#Check that we got what we asked for
	if not(mjd1 == dmjd[0] and mjd2 == dmjd[-1]):
		Debug('Bad data')
		sys.exit(0)
		
	# Save the data for post-mortems
	futcr = os.path.join(tmpDir,'utcr.{:d}.dat'.format(lastMJD))
	fout = open(futcr,'w')
	fout.write(resp.text)
	fout.close()
	Debug('Debug saved data as ' + futcr)
	
	# At this point, we can update
	fout = open(os.path.join(controlDir,lastUTCrFile),'w')
	fout.write('# MJD {:d} {:02d}:{:02d}:{:02d}\n'.format(mjdToday,dt.hour,dt.minute,dt.second))
	fout.write('MJD {:d}\n'.format(lastMJD))
	fout.close()
	
else: # TEMPORARY
	futcr = os.path.join(tmpDir,'utcr.{:d}.dat'.format(lastMJD))
	fin = open(futcr,'r')
	for l in fin:
		if re.match('\s*#',l): # ignore comments
			continue
		ldata = l.split()
		if len(ldata) == 2:
			try:
				dmjd.append(int(ldata[0]))
				dutck.append(float(ldata[1]))
			except:
				Log(logFile,'Bad data: ' + l)
	fin.close()

Debug('UTCr data: first = {:d}, last = {:d}'.format(dmjd[0],dmjd[-1]))

# Since we want do mucho maths, change to numpy arrays
mjd = np.array(dmjd)
utck = np.array(dutck)

# Check if steering is disabled eg for manual control
enableSteer = False # default is to NOT steer
if (os.path.exists(os.path.join(controlDir,enableSteerFile))):
	enableSteer = True
	
steerMsgs = ''

# TODO Check if there is an unprocessed steer

# This is temporary code for testing
FudgeData(mjd,utck)

if enableSteer or forceSteer:
	if enableSteer:
		Log(logFile,'Steering is ON')
	if forceSteer:
		Log(logFile,'Steering is FORCED')
		
	# Calculate the steering parameters
	
	
	# TODO Check when we last steered

	stopMJD  = mjd[-1]
	iStopMJD = len(mjd)-1
	startMJD = stopMJD
	iStartMJD = iStopMJD
	npts = 1

	for i in range(iStopMJD,0,-1): # search backwards
		if mjd[i] < stopMJD - 6:
			startMJD = dmjd[i+1]
			iStartMJD = i + 1
			npts = iStopMJD  - i
			break

	Debug('Fit: start = {:d}, stop = {:d}, npts = {:d}'.format(startMJD,stopMJD,npts))

	# Load the data masking file
	# This is just a list of MJDs
	mask = []
	fmask = os.path.join(controlDir,'mask.dat')
	if (os.path.exists(fmask)):
		Debug('Loading ' + fmask)
		fin = open(fmask,'r')
		for l in fin:
			if re.match('\s*#',l): # ignore comments
				continue
			matches = re.match(r'\s*(\d{5})\s*(\d{5})',l)
			if matches:
				Debug('Adding {} - {} to data mask'.format(matches.group(1),matches.group(2)))
				m1 = int(matches.group(1))
				m2 = int(matches.group(2))
				if (m1 > m2): # don't try to fix it - most likely it is a typo
					Debug('Bad mask: {:d} {:d}'.format(m1,m2)) 
				else:
					mask.append([m1,m2])
		fin.close
	
	# Calculate first MJD in the next UTCr report
	# It should be stopMJD + 1
	# nextStartMJD = 7*(ceil(stopMJD+1)/7);

	# Estimate the frequency offset to be applied to make the frequency offset == 0

	# If there are insufficient points for a fit then do not proceed
	# (Maybe we missed a UCTr reporting deadline)
	if npts < minFitPts:
		scheduleSteer = False
		msg = 'Insufficient points ({:d}) to estimate the mean frequency offset'.format(npts) 
		steerMsgs += msg + '<br>'
		Log(logFile,msg)

	# Filter data using the mask
	mjdFit = []
	utckFit = []
	nMasked = 0
	for im in range(iStartMJD,iStopMJD+1):
		masked = False
		for r in mask:
			if (mjd[im] >= r[0] and mjd[im] <= r[1]):
				masked = True
				nMasked += 1
				Debug('Masking {:d}'.format(mjd[im]))
		if not(masked): # so that we only add once per MJD
			mjdFit.append(mjd[im])
			utckFit.append(utck[im])
	Debug('Masked {:d} point(s)'.format(nMasked))
	
	# Steering algorithm, based on Chadsey et al 
	# A steer is constructed from two frequency adjustments
	# (1) Estimate the mean ffe using the last 7 days of UTCr and use this to zero the mean ffe
	# (2) Get the current (last day in UTCr) phase offset and apply a frequency offset to zero it over the next 7 days
	# The total steer is the weighted sum of the two
	
	coeff = np.polyfit( mjdFit - mjdFit[0],utckFit - utck[0],1) # units are ns/day
	meanfOffset = -coeff[0] # note the sign! Units are ns/day (the applied frequency steer will have opposite sign)
	freqStep = -meanfOffset
	
	Debug('Est. mean frequency offset: {:g} ns/day'.format(meanfOffset))
	Debug('->required frequency step:  {:g} ns/day'.format(freqStep))

	if abs(freqStep) > maxClockFOffset:
		scheduleSteer = False
		msg = 'Mean frequency offset exceeds limit (> {:g} ns/day)'.format(maxClockFOffset)
		steerMsgs += msg + '<br>'
		Log(logFile,msg)

	# Calculate the slew rate required to zero the offset
	currPhaseOffset = utckFit[-1];
	Debug('Current phase offset = {:g} ns'.format(currPhaseOffset))
	phaseSlew = currPhaseOffset/7;
	Debug('Est. phase offset slew: {:g} ns/day)'.format(phaseSlew))

	if abs(phaseSlew) > maxPhaseSlew:
		scheduleSteer = False
		msg = 'Required phase slew exceeds limit (> {:g} ns/day)'.format(maxPhaseSlew)
		steerMsgs += msg + '<br>'
		Log(logFile,msg)
	
	appliedOffset = freqGain * freqStep + phaseGain * phaseSlew
	# end of steering algorithm
	
else:
	Log(logFile,'Steering is OFF')

if (enableSteer and scheduleSteer) or forceSteer:
	steerFile = os.path.join(scheduledDir,dt.strftime('steer_%Y%m%d.dat'))
	try:
		ffe = appliedOffset*1.0E-9/86400.0
		fout = open(steerFile,'w')
		fout.write('# frequency steer {}\n'.format(dt.strftime('%Y%m%d')))
		fout.write('{:g}'.format(ffe))
		fout.close()
		Log(logFile,'Steer of {:g} scheduled'.format(ffe))
	except:
		Log(logFile,'Unable to create the steer file ' + steerFile)
		steerMsgs += '<strong>Unable to create the steer file ' + steerFile + '</strong>'
		# FIXME what to do next

# For the moment, we'll ignore the effect of gaps
statsMJD = [] # don't actually need this ATM
statsUTCK = []
mjdutck = zip(mjd,utck)

for m,u in mjdutck:
	masked = False
	for r in mask:
		if (m >= r[0] and m <= r[1]):
			masked = True
			break
	if not(masked):
		statsMJD.append(m)
		statsUTCK.append(u)

nData = len(statsMJD)
nTOTDEV = int(2.0*nData/3.0)
nTDEV   = int(nData/3)

# TOTDEV/...
t = np.arange(1,nTOTDEV,1)
(taus, devs, errors, ns) = allantools.totdev(statsUTCK,rate = 1.0,taus=t)

# TDEV/...
t = np.arange(1,nTDEV,1)
(tdtaus, tddevs, errors, ns) = allantools.tdev(statsUTCK,rate = 1.0,taus=t)

# Now make the report

# PHASE OFFSET
fig,ax = plt.subplots()
ax.plot(mjd,utck,label='all data')
ax.plot(statsMJD,statsUTCK,marker='.',label='stats data')
ax.set_ylabel('phase offset UTC (ns)')
ax.set_xlabel('MJD')
plt.legend()
ax.grid()

plotfile1 = os.path.join(tmpDir,'utcoffset.png')
fig.savefig(plotfile1)

if args.show:
	fig.show()

# TDEV
fig,ax = plt.subplots()

ax.loglog(tdtaus,tddevs,'.-')

ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
ax.set_ylabel('TDEV (ns)')
ax.set_xlabel('tau (days)')
ax.grid()

plotfile2 = os.path.join(tmpDir,'tdev.png')
fig.savefig(plotfile2)

if args.show:
	fig.show()

# TOTDEV
fig, ax = plt.subplots()
#ax3.set_title('UTC TOTDEV')
ax.loglog(taus,devs*1.0E-9/86400*1.0E14,'.-')

ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
ax.yaxis.set_major_formatter(mticker.ScalarFormatter())

ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax.set_ylabel('frac TOTDEV (10E-14 Hz/Hz)')
ax.set_xlabel('tau (days)')
ax.grid()

plotfile3 = os.path.join(tmpDir,'dev.png')
fig.savefig(plotfile3)

if args.show: 
	plt.show() # stops further execution

UTCstr = 'UTC(' + UTCID +')'
html =  '<html>'
html += '<head></head>'
html += '<body>'
html += '<H2>' + UTCstr + ' steering advisory for '+  dt.strftime('%Y-%m-%d') + '</H2>'
html += '<br>'

html += 'Current MJD is ' + str(mjdToday) + '<br>'
html += '<br>'

if not(enableSteer) and not(forceSteer):
	html += 'Automatic steering is currently disabled<br>'
	html += 'No steering parameters have been calculated<br>'
	html += '(re-enable by touching ~/control/' + enableSteerFile + ')<br>'
else:
	html += 'Est. mean frequency offset = {:g} ns/day<br>'.format(meanfOffset)
	html += 'Current phase offset = {:g} ns <br>'.format(currPhaseOffset)
	
	html += 'Est. frequency step  = {:g} ns/day <br>'.format(freqStep)
	html += 'Est. phase slew = {:g} ns/day <br>'.format(phaseSlew)
	
	html += 'Total (weighted) applied frequency offset = {:g} ns/day (ffe {:g})<br>'.format(appliedOffset,appliedOffset*1.0E-9/86400.0)
	if not(scheduleSteer):
		if forceSteer:
			html += '<div> <strong> STEERING WILL BE FORCED </strong> </div>'
		else:
			html += '<div> <strong> NO STEER WILL BE APPLIED </strong> </div>'
		html += steerMsgs 

html += '<H3>Time offset: UTCr - ' + UTCstr + '</H3>' 
html += '<img src="cid:offsetplot" alt="phase offset">'

html += '<H3>Time deviation</H3>' 
html += '<img src="cid:tdevplot" alt="time deviation">'

html += '<H3>Stability<H3>' 
html += '<img src="cid:devplot" alt="frequency deviation">'

html += '<br> <em>Generated by ' + appName + ' ' + VERSION + ' (' + user + '@' + socket.gethostname() + ' </em>)</br>\n' 
html += '</body>'
html += '</html>'

repBody = html
repBody = repBody.replace('cid:offsetplot',plotfile1) # hacky way to fix it up
repBody = repBody.replace('cid:tdevplot',plotfile2)
repBody = repBody.replace('cid:devplot',plotfile3)

repName = 'rep.{:d}{:02d}{:02d}'.format(yyyy,mm,dd)
repHTML = os.path.join(tmpDir,repName + '.html')
repPDF  = os.path.join(repDir,repName + '.pdf')
fout = open(repHTML,'w')
fout.write(repBody)
fout.close()

# Create a pdf from html
cmd = '/usr/bin/htmldoc --quiet --webpage -f ' + repPDF + ' ' + repHTML
os.system(cmd)

if email:
	sender = 'time@measurement.gov.au'

	msg = MIMEMultipart('related')
	msg['Subject'] = 'Weekly UTC steering advisory for ' + dt.strftime('%Y-%m-%d')
	msg['From'] = sender
	msg['To'] = recipients
	msg['Reply-To'] = sender

	body = MIMEText(html,'html')
	msg.attach(body)

	fp = open(plotfile1, 'rb')
	msgImage = MIMEImage(fp.read())
	fp.close()
	msgImage.add_header('Content-ID', '<offsetplot>')
	msg.attach(msgImage)

	fp = open(plotfile2, 'rb')
	msgImage = MIMEImage(fp.read())
	fp.close()
	msgImage.add_header('Content-ID', '<tdevplot>')
	msg.attach(msgImage)

	fp = open(plotfile3, 'rb')
	msgImage = MIMEImage(fp.read())
	fp.close()
	msgImage.add_header('Content-ID', '<devplot>')
	msg.attach(msgImage)

	# Send the message via local SMTP server.
	s = smtplib.SMTP('copperhead.in.measurement.gov.au')
	s.sendmail([recipients],['time@measurement.gov.au'], msg.as_string())
	s.quit()

# cleanup temporary files
if not(debug):
	os.unlink(plotfile1)
	os.unlink(plotfile2)
	os.unlink(plotfile3)
	#os.unlink(plotfile4)
	os.unlink(repHTML)
