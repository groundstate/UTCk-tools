#!/usr/bin/python3
#

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
# Check uploaded UTCr data
# This replaces UTCr_ftp_check.prl
#
# 
import argparse
import datetime
import hashlib
import math
import ftplib
import io
import numpy as np
import os
from pathlib import Path
import socket
import sys
import time

# This is where ottplib is installed
sys.path.append("/usr/local/lib/python3.6/site-packages")  # Ubuntu 18.04
sys.path.append("/usr/local/lib/python3.8/site-packages")  # Ubuntu 20.04
sys.path.append("/usr/local/lib/python3.10/site-packages") # Ubuntu 22.04
sys.path.append("/usr/local/lib/python3.12/site-packages") # Ubuntu 24.04

import ottplib as ottp

import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.image import MIMEImage

def md5sum(filePath):

	m = hashlib.md5()

	# Open the file in binary read mode ('rb')
	try:
		with open(filePath, 'rb') as f:
			# Read the file in chunks to be memory efficient for large files
			while True:
				data = f.read(8192) # Read in 8KB chunks
				if not data:
					break # Break if end of file is reached
				m.update(data)
		return m.hexdigest() # Return the hash as a hexadecimal string
	except Exception as e:
		return f"An error occurred: {e}"
			
VERSION = "0.4.1"
AUTHORS = "Michael Wouters"

NCHECK = 10        # number of MJDs to go back for check
HISTORY_LEN =  90  # number of days to plot

#
# Main 
#

appName = os.path.basename(sys.argv[0])
home =os.environ['HOME'] + '/'
user =os.environ['USER'] # remember to define this in the crontab
configFile = os.path.join(home,'etc/utcrcheck.conf')
tmpDir = os.path.join(home,'tmp')
recipients = 'time@measurement.gov.au'
sender = recipients
email = False
uploadDelay = 1 

tt = time.time()
mjdToday = int(tt/86400)+40587
dt = datetime.datetime.today()

parser = argparse.ArgumentParser(description='Check clock data uploaded to BIPM',
	formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--config','-c',help='use this configuration file',default=configFile)
parser.add_argument('--debug','-d',help='debug (to stderr)',action='store_true')
parser.add_argument('--show','-s',help='show plots',action='store_true')
parser.add_argument('--email','-e',help='send email',action='store_true')
parser.add_argument('--version','-v',action='version',version = os.path.basename(sys.argv[0])+ ' ' + VERSION + '\n' + 'Written by ' + AUTHORS)

args = parser.parse_args()

if (args.config):
	configFile = args.config
	if (not os.path.isfile(configFile)):
		ottp.ErrorExit(configFile + ' not found')

debug = args.debug
ottp.SetDebugging(debug)

email = args.email

if args.show:
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib.pyplot as plt
	import matplotlib.ticker as mticker
else:
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib as mplt # this (and the next line) stops warnings about being unable to connect to a display
	mplt.use('Agg') 
	import matplotlib.pyplot as plt
	import matplotlib.ticker as mticker
	
cfg=ottp.Initialise(configFile,['bipm:server','bipm:user','bipm:password',
	'data:file prefix','data:remote directory','data:local directory'])
server = cfg['bipm:server']
usr=cfg['bipm:user']
password=cfg['bipm:password']

remoteDir = cfg['data:remote directory']
localDir = ottp.MakeAbsolutePath(cfg['data:local directory'],home)
if 'data:tmp directory' in cfg:
	tmpDir = ottp.MakeAbsolutePath('data:tmp directory',home)	
prefix = cfg['data:file prefix']

if 'email:sender' in cfg:
	sender = cfg['email:sender']
if 'email:recipients' in cfg:
	recipients = cfg['email:sender']
if 'email:smtp server' in cfg:
	SMTPserver = cfg['email:smtp server']
	
# Useful for debugging problems
# ftp = ftplib.FTP(server)
# ftp.set_debuglevel(3)
# ftp.sanitize = lambda s: repr(s) # makes password visible
# ftp.login(user=usr, passwd=password)

# Construct the list of files to download and check
files = []
startMJD = mjdToday - uploadDelay - NCHECK + 1
stopMJD = mjdToday - uploadDelay
for m in range(startMJD, stopMJD + 1):
	dd = int(math.floor(m/1000))
	ddd = m - dd*1000
	fName = f'{prefix}__{dd}.{ddd:03d}' 
	files.append([fName,False,False,False])

checkOK = True
try:
	# Connect to the FTP server
	# The connection is in passive mode and with SSL off, bby default.
	with ftplib.FTP(server, usr, password) as ftp:
		ottp.Debug(f"Connected to {server}")

		# ftp.set_pasv(True) 
	
		ftp.cwd(remoteDir)
		
#		for i in range(0,len(files)):
		for i in range(0,len(files)):
			f = files[i]
			fp = Path(os.path.join(localDir,f[0]))
			ottp.Debug(f'Checking {f[0]}')
			if fp.is_file():
				ottp.Debug('Local: exists')
				files[i][1] = True
				localMD5sum = md5sum(fp)
				
				# To preserve line endings, download as binary adata
				bindata=bytearray()
				def callback(data): # appends downloaded data
					bindata.extend(data)
				try:      
					ftp.retrbinary(f'RETR {f[0]}',callback)
					files[i][2] = True
					ottp.Debug('Remote: exists')
					md5HashObject = hashlib.md5(bindata)
					remoteMD5sum = md5HashObject.hexdigest()
					files[i][3] = (remoteMD5sum == localMD5sum)
				except:
					files[i][2] = False
					ottp.Debug('Remote: missing')
			else:
				ottp.Debug(f'Local: missing')
				files[i][1] = False
				
		ftp.quit()
		
except ftplib.all_errors as e:
	print(f"FTP Error: {e}")
	checkOK = False

html =  '<html>\n'
html += '<head></head>\n'
html += '<body style="font-family: sans-serif">\n'
html += '<H2> UTCr data check for  '+  dt.strftime('%Y-%m-%d') + '</H2>\n'
html += '<br>\n'

html += 'The current MJD is ' + str(mjdToday) + '<br>\n'
html += '<br>\n'


if checkOK:
	html += '<H3>Check on uploaded clock data files:</H3>\n'
	for f in files:
		if not(f[1]):
			html += f' <div style="color: red;">{f[0]}    local file is missing</div>\n'
		elif not(f[2]):
			html += f' <div style="color: red;">{f[0]}    remote file is missing</div>\n'
		elif not(f[3]):
			html += f' <div style="color: red;">{f[0]}    local and renote files do not match</div>\n'
		else:
			html += f'<div>{f[0]}    OK<div>\n'
else:
	html += 'Check failed - network connectivity? <br>\n'

if checkOK: # make a plot of the local data
	# Use a hash table to store the data
	clkdata = {}
	clksteps = [] # save any detected clock steps that are within the reporting range

	for m in range(stopMJD - HISTORY_LEN,stopMJD + 1):
		dd = int(math.floor(m/1000))
		ddd = m - dd*1000
		fName = f'{prefix}__{dd}.{ddd:03d}'
		fp = Path(os.path.join(localDir,fName))
		# Format of each line is 
		# MJD LABID [CLKTYPE+CLKID  MEASUREMENT]
		if fp.is_file():
			fin = open(fp,'r')
			for l in fin:
				l = l.rstrip()
				# Test for a clock step line
				# This can be identified by there being an alphabetic character at position 41 (40 zero indexed)
				if len(l) >= 40 and l[40].isalpha(): # short circuited evaluation so safe
					# Save any clock steps that are recent
					mjd = int(l[0:5])
					if (mjd >= startMJD and mjd <= stopMJD):
						ottp.Debug(f'clock step:{l}')
						clksteps.append(l)
					continue
				d = l.split()
				if len(d)>= 4:
					mjd = d[0]
					for c in range(2,len(d),2):
						clkid = d[c][3:]
						if clkid[0] == '0':
							clkid = clkid[1:]
						if clkid in clkdata:
							clkdata[clkid][0].append(int(mjd))
							clkdata[clkid][1].append(float(d[c+1]))
						else:
							clkdata[clkid] = [[int(mjd)],[float(d[c+1])]]
			fin.close()

	fig,ax = plt.subplots()

	if args.show:
		fig.show()
	
	for clk in clkdata:
		tclk = np.array(clkdata[clk][0])
		uclk = np.array(clkdata[clk][1])
		coeff = np.polyfit(tclk,uclk,1) # units are ns/day
		foffset = -coeff[0] # note the sign!
		toffset =  coeff[1]
		ax.plot(tclk,uclk + foffset*tclk - toffset,label=clk,lw=1)
	ax.set_ylabel('UTC(AUS) - clock (ns)')
	ax.set_xlabel('MJD')
	plt.legend()
	ax.grid()

	plotfile = os.path.join(tmpDir,'utcrclocks.png')
	fig.savefig(plotfile)

	if args.show: 
		plt.show() # stops further execution
	
	# Add any clock step sto the report
	html += '<H3>Reported clock steps</H3>\n'
	if clksteps:
		for cs in clksteps:
			html += f'<div>{cs}</div>\n'	
	else:
		html += '<div>NONE</div>\n'
	html+= '<H3>Local clock data for the last 90 days</H3>\n'
	html += '<img src="cid:utcrclocks" alt="clock data plot">'
	html += '<div>Rate and mean offset have been removed in the plot</div>\n'
html += '<br> <em>Generated by ' + appName + ' ' + VERSION + ' (' + user + '@' + socket.gethostname() + ' </em>)</br>\n' 
html += '</body>\n'
html += '</html>\n'

ottp.Debug(html)

if email:
	
	msg = MIMEMultipart('related')
	msg['Subject'] = 'UTCr data check ' + dt.strftime('%Y-%m-%d')
	msg['From'] = sender
	msg['To'] = recipients
	msg['Reply-To'] = sender

	body = MIMEText(html,'html')
	msg.attach(body)

	fp = open(plotfile,'rb')
	msgImage = MIMEImage(fp.read())
	fp.close()
	msgImage.add_header('Content-ID', '<utcrclocks>')
	msg.attach(msgImage)
	
	# Send the message via local SMTP server.
	s = smtplib.SMTP(SMTPserver)
	s.sendmail(sender,recipients.split(','), msg.as_string())
	s.quit()
	
if not(debug):
	os.unlink(plotfile)
