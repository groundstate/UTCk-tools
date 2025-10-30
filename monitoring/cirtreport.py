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
import requests
import socket
import sys
import time

import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.image import MIMEImage

# This is where ottplib is installed
sys.path.append("/usr/local/lib/python3.6/site-packages")  # Ubuntu 18.04
sys.path.append("/usr/local/lib/python3.8/site-packages")  # Ubuntu 20.04
sys.path.append("/usr/local/lib/python3.10/site-packages") # Ubuntu 22.04
sys.path.append("/usr/local/lib/python3.12/site-packages") # Ubuntu 24.04

import ottplib as ottp

VERSION = "1.0.0"
AUTHORS = "Michael Wouters"

# ------------------------------------------
def Warn(msg):
	if (not args.nowarn):
		sys.stderr.write('WARNING! '+ msg+'\n')
	return
	
# ---------------------------------------------
def GetCircularT(lab,startMJD,stopMJD):
	ottp.Debug('Fetching Circular T data ...')
	try:
		r = requests.get(f'{httpRequest}scale=utc&lab={lab}&mjd1={startMJD}&mjd2={stopMJD}&outfile=txt')
	except:
		return None,None
	ottp.Debug('... done')
	lines = r.text.split('\r\n')
	mjds=[]
	utck = []
	for l in lines:
		l = l.strip()
		if not(l):
			continue
		if l[0] == '#':
			continue
		vals = l.strip().split()
		if (len(vals) == 3):
			mjds.append(int(vals[0]))
			utck.append(float(vals[1]))
	return mjds,utck
	
# ------------------------------------------

home =os.environ['HOME'] + '/'
user =os.environ['USER'] # remember to define this in the crontab
configFile = os.path.join(home,'etc/cirtreport.conf')
repDir  = os.path.join(home,'cirt/reports') 
historyLength = 12 # in months
tmpDir = os.path.join(home,'tmp')
recipients = 'time@measurement.gov.au'
email = False
lab = 'AUS'
httpRequest ='https://webtai.bipm.org/api/v1.0/get-data.html?'

tt = time.time()
mjdToday = int(tt/86400)+40587
dt = datetime.date.today()
yyyy = dt.year
mm   = dt.month
dd   = dt.day 

parser = argparse.ArgumentParser(description='Report on UTC(k) from Circular T data',
	formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--config','-c',help='use this configuration file',default=configFile)
parser.add_argument('--debug','-d',help='debug (to stderr)',action='store_true')
parser.add_argument('--show','-s',help='show plots',action='store_true')
parser.add_argument('--version','-v',action='version',version = os.path.basename(sys.argv[0])+ ' ' + VERSION + '\n' + 'Written by ' + AUTHORS)

args = parser.parse_args()

if (args.config):
	configFile = args.config
	if (not os.path.isfile(configFile)):
		ottp.ErrorExit(configFile + ' not found')

debug = args.debug
ottp.SetDebugging(debug)

cfg=ottp.Initialise(configFile,[])

if ('main:history' in cfg):
	historyLength = int(cfg['main:history'])

if ('main:lab' in cfg):
	lab = cfg['main:lab']

if ('main:email' in cfg):
	email = (cfg['main:email'].lower() == 'yes')
	
if ('main:email recipients' in cfg):
	recipients = cfg['main:email recipients']
	
if ('main:report path' in cfg):
	repDir = ottp.MakeAbsolutePath(cfg['main:report path'],home)

if ('main:tmp path' in cfg):
	tmpDir = ottp.MakeAbsolutePath(cfg['main:tmp path'],home) 
	
stopMJD  = mjdToday
startMJD = stopMJD - (historyLength+1)*31 # no need to be fussy here

ddmjd,ddutck = GetCircularT(lab,startMJD,stopMJD)
if not ddmjd:
	sys.exit('Failed to get CircularT')

startMJD = ddmjd[0]
stopMJD  = ddmjd[-1]

# Since we want to do maths change to numpy arrays
mjd = np.array(ddmjd)
utck = np.array(ddutck)

# Calculate the frequency offset
# Simple difference,  at midpoint
# Don't use np.diff(), in case there is missing data
fmjd = np.arange(mjd.size -1 ,dtype=float)
ffe  = np.arange(mjd.size -1 ,dtype=float)

for i in range(1,len(ddmjd)):
	fmjd[i-1] =(mjd[i] + mjd[i-1])/2.0
	ffe[i-1] = -1.0E-9*(utck[i] - utck[i-1]) / ((mjd[i] - mjd[i-1])*86400) # correct units Hz/Hz
	
minffe = np.min(ffe)
maxffe = np.max(ffe)

# Calculate the mean frequency from a linear fit to the data

coeff = np.polyfit(mjd - mjd[0],utck-utck[0],1) # units are ns/day
foffset = -coeff[0]*1.0E-9/86400 # note the sign!

# TOTDEV.
nd = int((stopMJD - startMJD)/10) # half the data
tau = 5*np.arange(1,nd)
(taus, devs, errors, ns) = allantools.totdev(utck,rate = 1.0/5.0,taus=tau)

# TDEV/...
nd = int((stopMJD - startMJD)/15)
tau = 5*np.arange(1,nd)			 
(tdtaus, tddevs, errors, ns) = allantools.tdev(utck,rate = 1.0/5.0,taus =tau)

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

ax.loglog(tdtaus,tddevs,'o-')

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
ax.loglog(taus,devs*1.0E-9/86400.0,'o-')
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

UTCstr = 'UTC(' + lab +')'
html =  '<html>'
html += '<head></head>'
html += '<body>'
html += '<H2>Report on '+ UTCstr + ' for ' + dt.strftime('%b %Y') + '</H2>'
html += '<br>'
# html += 'Most recent data is from Circular T issue ' + str(lastIssue) + '<br>' # note that lastIssue exists (but what if empty)
html += 'Last reported MJD is ' + str(ddmjd[-1]) + '<br>'
html += 'Current MJD is ' + str(mjdToday) + '<br>'
html += '<br>'
html += 'Last ' + str(ndays) + ' days <br>'
html += 'Average frequency offset = ' + '{:1.2E}'.format(foffset30d) + '<br>'
html += 'Frequency offset range (5d averages) ' + '[{:1.2E} {:1.2E}]'.format(minffe30,maxffe30) + '<br>'
html += '<br>'
html += f'Last {historyLength} months<br>'
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
	msg['Subject'] = f'Monthly UTC({lab}) report for ' + dt.strftime('%b %Y')
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

# Remove temporary files
if not(debug):
	os.unlink(plotfile1)
	os.unlink(plotfile2)
	os.unlink(plotfile3)
	os.unlink(plotfile4)
	os.unlink(repHTML)
