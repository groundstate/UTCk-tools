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
from fpdf import FPDF
import json
import numpy as np
import os
from pathlib import Path
import requests
import socket
import sqlite3
import sys

import smtplib
from email.message import EmailMessage

try:
	import ottplib as ottp
except ImportError:
	sys.exit('ERROR: Must install ottplib\n eg openttp/software/system/installsys.py -i ottplib')

VERSION = '0.6.0'
AUTHORS = 'Michael Wouters'

# ------------------------------------------
def ShowVersion():
	print (os.path.basename(sys.argv[0])+" "+VERSION)
	print ('Written by ' + AUTHORS)
	return

# ---------------------------------------------
def FetchBIPMBUTC(gnss,startMJD,stopMJD):
	
	ottp.Debug(f'Fetching Circular bUTC prediction data  for {gnss} from BIPM');
	if debug:
		#TEMPORARY CODE
		with open(os.path.join(root,'database/butc_gps.json'), 'r', encoding='utf-8') as file:
			d = json.load(file)
	else:
		req = f'{httpRequest}scale=b_{gnss}&mjd1={startMJD}&mjd2={stopMJD}&outfile=json'
		ottp.Debug(req)
		try:
			r = requests.get(req,verify=rootCert)
		except:
			ottp.Debug('http request failed')
			return None,None,None
		
		d = json.loads(r.text)
		if not(d['errorcode'] == 0):
			ottp.Debug(f"http request return error {d['errorcode']}")
			return None,None,None
		
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
user =os.environ['USER']
host = socket.gethostname()
root = home

httpRequest = 'https://webtai.bipm.org/api/v1.0/get-data.html?'
rootCert = None # if you use an empty string, this will skip SSL verification which is a bad thing
SMTPserver = 'localhost'
lab = 'k'

checkingPath = os.path.join(root,'report/checking')

configFile = os.path.join(root,'etc/butc.conf')

if ottp.LibMajorVersion() >= 0 and ottp.LibMinorVersion() < 2: 
	sys.exit('Need ottplib minor version >= 2')

examples='TO DO'
parser = argparse.ArgumentParser(description='Provides information for checking a monthly report',
	formatter_class=argparse.RawDescriptionHelpFormatter,epilog=examples)
parser.add_argument('--config','-c',help='use an alternate configuration file',default=configFile)

parser.add_argument('--debug','-d',help='debug (to stderr)',action='store_true')
parser.add_argument('--display',help='display plots',action='store_true')
parser.add_argument('--email',help='email report for checking',action='store_true')
parser.add_argument('--force',help='force creation of a new report',action='store_true')
parser.add_argument('--version','-v',help='show version and exit',action='store_true')
parser.add_argument('--month','-m',help='create a report for the given month 1..12')
parser.add_argument('--year','-y',help='create a report for the given year')

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

cfg=ottp.Initialise(configFile,['main:gnss','report:checking path','email:recipients','email:sender'])

gnss = cfg['main:gnss'].split(',')
gnss = [g.strip() for g in gnss] 

if 'main:lab' in cfg:
	lab = cfg['main:lab']
	
if 'paths:root' in cfg:
	root = cfg['paths:root']
	# If the root is not absolute, prepend the user's home directory
	tmpPath = Path(root)
	if not tmpPath.is_absolute():
		root = os.path.join(home,root)

checkingPath = ottp.MakeAbsolutePath(cfg['report:checking path'],root)

recipients =  cfg['email:recipients']
sender = cfg['email:sender']

if ('email:smtp server' in cfg):
	SMTPserver = cfg['email:smtp server']

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

if (args.month is not None) != (args.year is not None):
	ottp.ErrorExit('You need to specify --month and --year')
elif (args.month and args.year):
	yyyy = int(args.year)
	mm = int(args.month)
	
	mmStart = mm
	yyyyStart = yyyy
	mmStop = mm
	yyyyStop = yyyy
	
	_,lastDayOfMonth = calendar.monthrange(yyyyStop,mmStop)
	
# If the reporting period has been set on the command line, then we 
# check if anything needs to be done
if not(args.month):
	
	# Do we need to do this ?
	# Is last month's PDF plot already there?
	fPath  = Path(os.path.join(checkingPath,f'butcplot_{yyyyStop}{mmStop:02d}.pdf'))
	if fPath.is_file() and (not args.force):
		ottp.Debug(f'butcplot_{yyyyStop}{mmStop:02d}.pdf already created - nothing to do')
		sys.exit(0)
		
	# Are new UTC data available for release?
	startMJD = ottp.MJD(datetime(yyyyStart,mmStart,1,0,0,0,tzinfo=timezone.utc).timestamp())
	stopMJD =  ottp.MJD(datetime(yyyyStop,mmStop,lastDayOfMonth,0,0,0,tzinfo=timezone.utc).timestamp())

	# Search for unreleased data
	ottp.Debug(f'Connecting to the database {db}')
	dbc = sqlite3.connect(db) # this opens a connection to the database
	curs = dbc.cursor()  
	firstUnreleased = -1
	lastUnreleased= -1
	for m in range (startMJD, stopMJD+1):
		r = curs.execute(f'SELECT release_utc FROM butcgnss WHERE mjd={m};')
		x = r.fetchone()
		if x:
			if x[0] == 0:
				lastUnreleased = m
				if firstUnreleased < 0:
					firstUnreleased = m
	ottp.Debug(f'Unreleased data {firstUnreleased}-{lastUnreleased}')

	curs.close()
	dbc.close()
	
	if lastUnreleased == -1 and not(args.force):
		ottp.Debug('No unreleased data')
		sys.exit(0)

ottp.Debug(f'Month range {mmStart:02d}-{yyyyStart}->{mmStop:02d}-{yyyyStop}')

startMJD = ottp.MJD(datetime(yyyyStart,mmStart,1,0,0,0,tzinfo=timezone.utc).timestamp())
stopMJD =  ottp.MJD(datetime(yyyyStop,mmStop,lastDayOfMonth,0,0,0,tzinfo=timezone.utc).timestamp())

ottp.Debug(f'Processing for {startMJD} to {stopMJD}')

# Make reports for each GNSS

# Make a plot for each GNSS - do this first so they can be interleaved with the data
fig, axs = plt.subplots(1,1,figsize=[8,8],squeeze=False ) #unsqueeze so we always get an array of axes
#title = 'bUTC_GNSS check for ' + calendar.month_name[mmStop]+ f' {yyyyStop}\n'
#title += os.path.basename(sys.argv[0])+ ' v' + VERSION   + ' run ' + dt.strftime('%Y-%m-%d %H:%M:%S') + '\n'
#fig.suptitle(title,ha='left',x=0.1)

plotIndex = 0

plots={}
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

	plotFile = os.path.join(checkingPath,f'{g}_{yyyyStop}{mmStop:02d}.png')
	plt.savefig(plotFile,format='png',) # do this before show()
	plots[g]=plotFile

dbc = sqlite3.connect(db) # this opens a connection to the database
curs = dbc.cursor()   

missingData = '*' 

pdf = FPDF()
pdf.set_font('helvetica', size=10)
pdf.set_margins(20, 20, 20) # in mm

pdf.add_page()
html = f'<H2> bUTC_GNSS checking report for {calendar.month_name[mmStop]} {yyyyStop} </H2>'
html += '<br> Generated by </br>'
html += f'<br> {user}@{host} {os.path.basename(sys.argv[0])} {VERSION}<br>'
pdf.write_html(html)

for gi in range(0,len(gnss)):
	g = gnss[gi]
	pdf.add_page()
	
	html = f'<h2 align="center">{g}</h2>'
	html += '<table>'
	html += '<tr>'
	html += f'<th>MJD</th> <th>UTC({lab})-bUTC</th> <th>U</th> <th>UTC-bUTC</th> <th>U</th>'
	html += '</tr>'
	for m in range(firstUnreleased,lastUnreleased+1):
		
		html += f'<tr align="right"><td align="center">{m:5d}</td>' # note that inline CSS is not supported
		
		r = curs.execute(f'SELECT UTCk_{g},UTCk_{g}_u,UTC_{g},UTC_{g}_u,release_utc from butcgnss where MJD={m};')
		x = r.fetchone()
		
		if x:
			if x[0]==None: # no data
				html += f'<td>{missingData}</td> <td>{missingData}</td> <td>{missingData}</td> <td>{missingData}</td>'
			elif x[2] == None : # no UTC data
				outputLine += '{:>9.2f} {:>7.2f} {:>9} {:>7}'.format(x[0],x[1],missingData,missingData)
			else: # Yay got it all
				#outputLine += '{:>9.2f} {:>7.2f} {:>9.2f} {:>7.2f}'.format(x[0],x[1],x[2],x[3])
				html += f'<td>{x[0]:>9.2f}</td> <td>{x[1]:>7.2f}</td> <td>{x[2]:>9.2f}</td> <td>{x[3]:>7.2f}</td>'
		else:
			outputLine += f'{missingData:>9} {missingData:>7} {missingData:>7} {missingData:>9}'
		
		html += '</tr>'
		
	html += '</table>'
	pdf.write_html(html)
	#pdf.add_page()
	pdf.image(plots[g],w=160)
		
curs.close()
dbc.close()

pdfFile = os.path.join(checkingPath,f'butc_{yyyyStop}{mmStop:02d}.pdf')
pdf.output(pdfFile)

if args.email:
	ottp.Debug('Sending email')
	msg =  EmailMessage()
	msg['Subject'] = 'bUTC_GNSS checking data for ' + calendar.month_name[mmStop]+ f' {yyyyStop}'
	msg['From'] = sender
	msg['To'] = recipients
	msg['Reply-To'] = sender

	msg.set_content('Save the attachments ...')
	
	# Add attachments
	with open(pdfFile,'rb') as fin:
		fData = fin.read()
		msg.add_attachment(fData,maintype='application',subtype='pdf',filename=Path(pdfFile).name)
	
	# Send the message via local SMTP server.
	s = smtplib.SMTP(SMTPserver)
	s.sendmail(sender,recipients.split(','), msg.as_string())
	s.quit()
	
if args.display:
	plt.show()

if not debug:
	for g in gnss:
		os.unlink(plots[g])


