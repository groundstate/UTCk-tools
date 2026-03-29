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
from pathlib import Path
import requests
import sqlite3
import sys

import smtplib
from email.message import EmailMessage

try:
	import ottplib as ottp
except ImportError:
	sys.exit('ERROR: Must install ottplib\n eg openttp/software/system/installsys.py -i ottplib')

VERSION = '0.2.0'
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
	with open(os.path.join(home,'database/butc_gps.json'), 'r', encoding='utf-8') as file:
		d = json.load(file)

	return np.array(d['data'][0]['x']),np.array(d['data'][0]['y']),np.array(d['data'][0]['unc'],dtype=float)

def ReadLabBUTC(db,g,startMJD,stopMJD):
	m  = []
	d  = []
	u  = []
	ottp.Debug(f'Connecting to the database {db}')
	dbc = sqlite3.connect(db) # this opens a connection to the database
	curs = dbc.cursor()   

	for mjd in range(startMJD,stopMJD+1):
		r = curs.execute(f'SELECT UTC_{g},UTC_{g}_u,release_utc from butcgnss where MJD={mjd};')
		x = r.fetchone()
		if x:
			if not (x[0] == None):
				m.append(mjd)
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

emailRecipients = 'time@measurement.gov.au'
emailSender = 'time@measurement.gov.au'
SMTPserver = 'localhost'

checkingPath = os.path.join(root,'report/checking')

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
parser.add_argument('--email',help='email data',action='store_true')
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

cfg=ottp.Initialise(configFile,['main:gnss','checking:path'])

gnss = cfg['main:gnss'].split(',')
gnss = [g.strip() for g in gnss] 

checkingPath = ottp.MakeAbsolutePath(cfg['checking:path'],root)

if ('main:email recipients' in cfg):
	recipients = cfg['main:email recipients']

if ('main:email sender' in cfg):
	emailSender = cfg['main:email sender']

if ('main:smtp server' in cfg):
	SMTPserver = cfg['main:smtp server']

if ('main:root certificate' in cfg):
	rootCert= cfg['main:root certificate']
	
if ('database:file' in cfg):
	db = ottp.MakeAbsoluteFilePath(cfg['database:file'],root,os.path.join(home,'database'))

dt   = datetime.now(tz=timezone.utc) # get date in UTC
mm   = dt.month
yyyy = dt.year
dd   = dt.day

mmStart = mm - 3 # need to check the previous two months
yyyyStart = yyyy
if mmStart <= 0:
	mmStart += 12
	yyyyStart = yyyy - 1
	
mmStop = mm - 1 
yyyyStop = yyyy
if mmStop <= 0:
	mmStop += 12
	yyyyStop = yyyy - 1
	
_,lastDayOfMonth = calendar.monthrange(yyyyStop,mmStop)

# If this has not been run manually then we do checks
if not(args.month):
	
	# Do we need to do this ?
	# Is last month's PDF plot already there?
	fPath  = Path(os.path.join(checkingPath,f'butcplot_{yyyyStop}{mmStop:02d}.pdf'))
	if fPath.is_file():
		ottp.Debug(f'butcplot_{yyyyStop}{mmStop:02d}.pdf already created - nothing to do')
		sys.exit(0)
		
	# Are new UTC data available for release?
	startMJD = ottp.MJD(datetime(yyyyStop,mmStop,1,0,0,0,tzinfo=timezone.utc).timestamp())
	stopMJD =  ottp.MJD(datetime(yyyyStop,mmStop,lastDayOfMonth,0,0,0,tzinfo=timezone.utc).timestamp())

	# Search backwards through the past month for unreleased data
	ottp.Debug(f'Connecting to the database {db}')
	dbc = sqlite3.connect(db) # this opens a connection to the database
	curs = dbc.cursor()  

	latestNewData = -1
	for m in range(stopMJD,startMJD-1,-1):
		r = curs.execute(f'SELECT release_utc FROM butcgnss WHERE mjd={m};')
		x = r.fetchone()
		if x:
			if x[0] == 0:
				ottp.Debug(f'Most recent unreleased data is for {m}')
				latestNewData = m
				break
	curs.close()
	dbc.close()
	
	if latestNewData == -1:
		ottp.Debug('No unreleased data')
		sys.exit(0)
	
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
	_,lastDayOfMonth = calendar.monthrange(yyyyStop,mmStop)
	
ottp.Debug(f'{mmStart} {yyyyStart},{mmStop} {yyyyStop}')

startMJD = ottp.MJD(datetime(yyyyStart,mmStart,1,0,0,0,tzinfo=timezone.utc).timestamp())
stopMJD =  ottp.MJD(datetime(yyyyStop,mmStop,lastDayOfMonth,0,0,0,tzinfo=timezone.utc).timestamp())

ottp.Debug(f'Processing for {startMJD} to {stopMJD}')

# Make a plot for each GNSS
# One page should be fine
fig, axs = plt.subplots(len(gnss),1,figsize=[8,12],squeeze=False) #unsqueeze so we always get an array of axes
title = 'bUTC_GNSS check for ' + calendar.month_name[mmStop]+ f' {yyyyStop}\n'
title += os.path.basename(sys.argv[0])+ ' v' + VERSION   + ' run ' + dt.strftime('%Y-%m-%d %H:%M:%S') + '\n'
fig.suptitle(title,ha='left',x=0.1)

plotIndex = 0

for g in gnss:
	m,d,u = FetchBIPMBUTC(g,startMJD,stopMJD) 
	
	axs[plotIndex,0].fill_between(m,d-u,d+u,color='palegreen')
	axs[plotIndex,0].plot(m,d,color='g')
	axs[plotIndex,0].scatter(m,d,s=9,marker='o',color='g',label='Circular T')
	
	m,d,u = ReadLabBUTC(db,g,startMJD,stopMJD)
	
	axs[plotIndex,0].fill_between(m,d-u,d+u,color='lightskyblue',alpha=0.5)
	axs[plotIndex,0].plot(m,d,color='b')
	axs[plotIndex,0].scatter(m,d,s=9,marker='o',color='b',label = 'local')

	axs[plotIndex,0].set_ylabel('ns')
	axs[plotIndex,0].set_xlabel('MJD')
	axs[plotIndex,0].set_title(f'UTC - bUTC_{g}')
	axs[plotIndex,0].grid()
	
	axs[plotIndex,0].legend()
	
	plotIndex += 1

plotFile = os.path.join(checkingPath,f'butcplot_{yyyyStop}{mmStop:02d}.pdf')
plt.savefig(plotFile,format='pdf',) # do this before show()

if args.email:
	
	msg =  EmailMessage()
	msg['Subject'] = 'bUTC_GNSS checking data for ' + calendar.month_name[mmStop]+ f' {yyyyStop}'
	msg['From'] = emailSender
	msg['To'] = emailRecipients
	msg['Reply-To'] = emailSender

	msg.set_content('Save the attachments ...')
	
	# Add attachments
	with open(plotFile,'rb') as fin:
		fData = fin.read()
		msg.add_attachment(fData,maintype='application',subtype='pdf',filename=Path(plotFile).name)
	
	# Pick up the monthly reports
	# TODO
	
	# Send the message via local SMTP server.
	s = smtplib.SMTP(SMTPserver)
	s.sendmail(emailSender,emailRecipients.split(','), msg.as_string())
	s.quit()
	
if args.display:
	plt.show()
	


