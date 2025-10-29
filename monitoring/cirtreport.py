#!/usr/bin/python3
#

#
# The MIT License (MIT)
#
# Copyright (c) 2018 Michael J. Wouters
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
# Report on UTC(AUS) using Circular T
# 

import allantools
import argparse
import datetime
import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mticker
import os
import re
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

VERSION = "0.1.0"
AUTHORS = "Michael Wouters"

# ------------------------------------------
def Debug(msg):
	if (debug):
		sys.stderr.write(msg+'\n')
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
		if (not cfg.has_key(k)):
			ErrorExit("The required configuration entry " + k + " is undefined")
		
	return cfg

# ------------------------------------------
def ReadCIRT(fname):
	# File f is presumed to be readable
	# Read header until we get the MJD line
	Debug('Opening ' + fname)
	fin = open(fname,'r')
	mjds=[]
	utck = []
	for l in fin:
		if (re.search(r'^\s*MJD',l)): # lazy
			args = l.strip().split()
			for i in range(1,len(args)-3):
				mjds.append(int(args[i]))
			break
	
	lab = r'^' + UTCID
	for l in fin:
		if (re.search(lab,l)):
			# Remove lab name and location for ease of parsing
			indx = l.find(')')
			ll = l[indx+1:]
			# Occasionally there are notes at the end
			indx = ll.find('(')
			if indx:
				ll = ll[0:indx]
			args = ll.strip().split()
			for i in range(0,len(args)-3):
				if (args[i] == '-'): # missing data
					utck.append(None) # tag for cleanup
				else:
					utck.append(float(args[i]))
			break
	
	fin.close()
	
	if (len(mjds) == len(utck)):
		# Remove data tagged with 'None'
		for i in range(0,len(mjds)):
			if utck[i] == None:
				mjds[i] = None
		mjds = [x for x in mjds if x is not None]
		utck = [x for x in utck if x is not None]
		return [mjds,utck]
	else:
		Debug('Mismatch in ' + fname)
		return [[],[]]
	
# ------------------------------------------

debug = False

UTCID = 'AUS'

home =os.environ['HOME'] + '/'
user =os.environ['USER'] # remember to define this in the crontab
configFile = os.path.join(home,'etc/cirtreport.conf')
cirtDir = os.path.join(home,'cirt') + os.sep
repDir  = os.path.join(home,'cirt/reports') 
historyLength = 6
tmpDir = os.path.join(home,'tmp')
recipients = 'time@measurement.gov.au'
email = True

tt = time.time()
mjdToday = int(tt/86400)+40587
dt = datetime.date.today()
yyyy = dt.year
mm   = dt.month - 1
if (mm == 0):
	mm = 12
	yyyy = yyyy - 1
dt=datetime.date(yyyy,mm,1)


parser = argparse.ArgumentParser(description='Report on UTC(k) from Circular T data',
	formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--config','-c',help='use this configuration file',default=configFile)
parser.add_argument('--debug','-d',help='debug (to stderr)',action='store_true')
parser.add_argument('--show','-s',help='show plots',action='store_true')
parser.add_argument('--version','-v',action='version',version = os.path.basename(sys.argv[0])+ ' ' + VERSION + '\n' + 'Written by ' + AUTHORS)

args = parser.parse_args()

debug = args.debug
configFile = args.config

cfg = Initialise(configFile)

if ('main:history' in cfg):
	historyLength = int(cfg['main:history'])

if ('main:utc id' in cfg):
	UTCID = cfg['main:utc id']

if ('main:email recipients' in cfg):
	recipients = cfg['main:email recipients']
	
if ('main:reports' in cfg):
	repDir = os.path.join(home,cfg['main:reports']) 
	
# Find the first and last Circular T files we have, by issue
ct=glob.glob(cirtDir + 'cirt.*')
if (0==len(ct)):
	ErrorExit('No Circular T files')
match = re.search('cirt\.(\d+)$',ct[0])
firstIssue= int(match.group(1))
lastIssue = int(match.group(1))

for c in ct:
	match = re.search('cirt\.(\d+)$',c)
	issue = int(match.group(1))
	if (issue < firstIssue):
		firstIssue = issue
	if (issue > lastIssue):
		lastIssue = issue

stop = lastIssue
start = lastIssue - historyLength
if (start < firstIssue):
	start = firstIssue
mjds = []
utck = []

for issue in range(start,stop+1):
	fname = cirtDir + 'cirt.' + str(issue)
	if (not os.path.isfile(fname)):
		Warn(fname + ' is missing')
		continue
	
	(tmjd,tutck) = ReadCIRT(fname)
	mjds = mjds + tmjd
	utck = utck + tutck

if len(mjds) ==0 or len(utck)==0:
	ErrorExit('No data')
	
# Now eliminate duplicates, retaining the most recent
ddmjd = [mjds[0]]
ddutck= [utck[0]]
for i in range(1,len(mjds)):
	if (mjds[i] == mjds[i-1]):
		# overwrite the last UTC offset 
		ddutck[-1] = utck[i]
	else: # append new values
		ddmjd = ddmjd + [mjds[i]]
		ddutck = ddutck + [utck[i]]

# Since we want do maths change to numpy arrays
mjd = np.array(ddmjd)
utck = np.array(ddutck)

# Calculate the frequency offset
# Simple difference,  at midpoint
# Don't use np.diff(), in case there is missing data
fmjd = np.arange(mjd.size -1 ,dtype=float)
ffe = np.arange(mjd.size -1 ,dtype=float)

for i in range(1,len(ddmjd)):
	fmjd[i-1] =(mjd[i] + mjd[i-1])/2.0
	ffe[i-1] = -1.0E-9*(utck[i] - utck[i-1]) / ((mjd[i] - mjd[i-1])*86400) # correct units Hz/Hz
	
minffe = np.min(ffe)
maxffe = np.max(ffe)

# Calculate the mean frequency from a linear fit to the data

coeff = np.polyfit(mjd - mjd[0],utck-utck[0],1) # units are ns/day
foffset = -coeff[0]*1.0E-9/86400 # note the sign!

# ADEV/...
(taus, devs, errors, ns) = allantools.totdev(utck,rate = 1.0/(5.0*86400.0))

# TDEV/...
(tdtaus, tddevs, errors, ns) = allantools.tdev(utck,rate = 1.0/(5.0*86400.0))

# Now make the report

#f,(ax1,ax2,ax3)= plt.subplots(3,sharex=False,figsize=(8,11))
#f.suptitle('UTC - ' + UTClab)

# TIME OFFSET

fig,ax = plt.subplots()
ax.plot(mjd,utck,marker='.')
mjdfit = np.array([0,mjd[-1]-mjd[0]])
utckfit = np.polyval(coeff,mjdfit)
                 
ax.plot(mjdfit + mjd[0],utckfit+utck[0],linestyle='dashed')
#ax1.set_title('UTC time offset')
ax.set_ylabel('time offset (ns)')
ax.set_xlabel('MJD')
ax.grid()

plotfile1 = os.path.join(tmpDir,'utcoffset.png')
fig.savefig(plotfile1)

if args.show:
	fig.show()

# TDEV

fig,ax = plt.subplots()

ax.loglog(tdtaus/86400,tddevs,'o-')

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

# FFE

fig,ax = plt.subplots()
ax.plot(fmjd,ffe*1.0E14,marker='.')
ax.plot([fmjd[0],fmjd[-1]],[foffset*1.0E14,foffset*1.0E14],linestyle='dashed')
#ax2.set_title('UTC frequency offset')
ax.set_ylabel(r'ffe ($10^{-14}$ Hz/Hz)')
ax.set_xlabel('MJD')
ax.grid()

plotfile3 = os.path.join(tmpDir,'ffe.png')
fig.savefig(plotfile3)

if args.show:   
	fig.show()

# TOTDEV

fig, ax = plt.subplots()
#ax3.set_title('UTC TOTDEV')
ax.loglog(taus/86400,devs*1.0E-9,'o-')
ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax.set_ylabel('frac TOTDEV (Hz/Hz)')
ax.set_xlabel('tau (days)')
ax.grid()

plotfile4 = os.path.join(tmpDir,'dev.png')
fig.savefig(plotfile4)

if args.show: 
	plt.show() # stops further execution

# Stats for the last 30 days
mjd30 = mjd[-7:]
utck30 = utck[-7:]
coeff30d = np.polyfit(mjd30 - mjd30[0],utck30 - utck30[0],1) # units are ns/day
foffset30d = -coeff30d[0]*1.0E-9/86400 # note the sign!
minffe30 = np.min(ffe[-7:])
maxffe30 = np.max(ffe[-7:])
ndays = mjd30[-1] - mjd30[0]

UTCstr = 'UTC(' + UTCID +')'
html =  '<html>'
html += '<head></head>'
html += '<body>'
html += '<H2>Report on '+ UTCstr + ' for ' + dt.strftime('%b %Y') + '</H2>'
html += '<br>'
html += 'Most recent data is from Circular T issue ' + str(lastIssue) + '<br>' # note that lastIssue exists (but what if empty)
html += 'Last reported MJD is ' + str(ddmjd[-1]) + '<br>'
html += 'Current MJD is ' + str(mjdToday) + '<br>'
html += '<br>'
html += 'Last ' + str(ndays) + ' days <br>'
html += 'Average frequency offset = ' + '{:1.2E}'.format(foffset30d) + '<br>'
html += 'Frequency offset range (5d averages) ' + '[{:1.2E} {:1.2E}]'.format(minffe30,maxffe30) + '<br>'
html += '<br>'
html += 'Last 12 months<br>'
html += 'Average frequency offset = ' + '{:1.2E}'.format(foffset) + '<br>'
html += 'Frequency offset range (5d averages) ' + '[{:1.2E} {:1.2E}]'.format(minffe,maxffe) + '<br>'

html += '<H3>Time offset: UTC - ' + UTCstr + '</H3>' 
html += '<img src="cid:offsetplot" alt="time offset">'

html += '<H3>Time deviation</H3>' 
html += '<img src="cid:tdevplot" alt="time deviation">'

html += '<H3>Frequency offset</H3>' 
html += '<img src="cid:ffeplot" alt="frequency offset">'

html += '<H3>Stability<H3>' 
html += '<img src="cid:devplot" alt="frequency deviation">'

html += '<br> <em>Generated by cirtreport.py ' + VERSION + ' (' + user + '@' + socket.gethostname() + ' </em>)</br>\n' 
html += '</body>'
html += '</html>'

repBody = html
repBody = repBody.replace('cid:offsetplot',plotfile1) # hacky way to fix it up
repBody = repBody.replace('cid:tdevplot',plotfile2)
repBody = repBody.replace('cid:ffeplot',plotfile3)
repBody = repBody.replace('cid:devplot',plotfile4)

repName = 'rep.{:d}{:02d}'.format(yyyy,mm)
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
	msg['Subject'] = 'Monthly UTC(AUS) report for ' + dt.strftime('%b %Y')
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
	msgImage.add_header('Content-ID', '<ffeplot>')
	msg.attach(msgImage)

	fp = open(plotfile4, 'rb')
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
	os.unlink(plotfile4)
	os.unlink(repHTML)
