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
# Executes a scheduled steer
# 

import argparse
import datetime
import glob
import os
import re
import serial
import shutil
import subprocess
import sys
import time

import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.image import MIMEImage

# This is where ottplib is installed
sys.path.append("/usr/local/lib/python3.6/site-packages") # Ubuntu 18.04
sys.path.append("/usr/local/lib/python3.8/site-packages") # Ubuntu 20.04
sys.path.append("/usr/local/lib/python3.10/site-packages") # Ubuntu 22.04

import ottplib

VERSION = "0.0.6"
AUTHORS = "Michael Wouters"

BaudRates = [9600,14400,19200,28800,38400,57600,115200]

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
		flog.write('{} {}:{}\n'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.gmtime()),appName,msg))
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
	reqd = ['main:lock file']
	for k in reqd:
		if (not k in cfg):
			ErrorExit("The required configuration entry " + k + " is undefined")
		
	return cfg

#-----------------------------------------------------------------------------
def Cleanup():
	# Hmm ugly globals
	ottplib.RemoveProcessLock(lockFile)
	if (not serport==None):
		serport.close()
		subprocess.check_output(['/usr/local/bin/lockport','-r',port])

#-----------------------------------------------------------------------------
def HROGCmd(serport,cmd):
	Debug('Sending ' + cmd)
	ba = bytearray() # stupid unicode
	cmd = cmd + '\r'
	ba.extend(cmd.encode())
	serport.write(ba)
	
#-----------------------------------------------------------------------------
def HROGQuery(serport,cmd):
	
	for i in range(0,5): # have a few goes
		HROGCmd(serport,cmd)
		d = serport.readline().decode('utf-8').strip()  # assuming  here that ECHO is off
		if re.search(cmd,d):
			Debug('Received ' + d)
			return d.split() # return a list of strings since usual response is CMD VALUE UNIT 
	return []

#-----------------------------------------------------------------------------
def GetHROGSettings(serport):
	stoffs = None
	ffof = None
	phase = None
	ret = HROGQuery(serport,'TOFFS?') # time offset TOFFS?
	# Don't convert return strings so we always report them as the HROG reports them
	if ret:
		stoffs = ret[1] # WARNING assuming units are always ns
	ret = HROGQuery(serport,'FFOF?') # current fractional frequency offset FFOF?
	if ret:
		ffof = ret[1]
	ret = HROGQuery(serport,'PHAS?')# phase offset PHAS?
	if ret:
		 phas = ret[1]
	return [stoffs,ffof,phas]

# ------------------------------------------
# Main

appName = os.path.basename(sys.argv[0])
debug = False
email = True
recipients = 'Michael.Wouters@measurement.gov.au'
UTCID = 'AUS'

# Limits on the applied frequency step
maxFSTEP = 1.0E-13 # a bit less than our accredited uncertainty
maxFFOF  = 5.0E-13 # maximum tolerated fractional frequency offset

tt = time.time()
mjdToday = int(tt/86400)+40587
Debug('MJD today is {:d}'.format(mjdToday))
dt = datetime.datetime.today()

home =os.environ['HOME'] + '/'
user =os.environ['USER'] # remember to define this in the user's crontab
configFile = os.path.join(home,'etc/utcsteer.conf')
ussConfigFile = os.path.join(home,'etc/utcsteersched.conf')
logDir = os.path.join(home,'logs')
controlDir = os.path.join(home,'control')
scheduledDir = os.path.join(controlDir,'scheduled_steer')
processedDir = os.path.join(controlDir,'processed_steer')

SMTPserver = 'copperhead.in.measurement.gov.au'
emailSender = 'time@measurement.gov.au'

parser = argparse.ArgumentParser(description='Executes a steer calculated and scheduled by utcsteersched.py',
	formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--config','-c',help='use an alternate configuration file',default=configFile)
parser.add_argument('--debug','-d',help='debug',action='store_true')
parser.add_argument('--dryrun',help='do everything except steer and move the steer file',action='store_true')
parser.add_argument('--show','-s',help='show HROG settings',action='store_true')
parser.add_argument('--version','-v',action='version',version = os.path.basename(sys.argv[0])+ ' ' + VERSION + '\n' + 'Written by ' + AUTHORS)

args = parser.parse_args()

debug = args.debug

configFile = args.config;

if (not os.path.isfile(configFile)):
	ErrorExit(configFile + ' not found')

cfg=Initialise(configFile)

# Get some settings from utcsteersched.conf
ussCfg = ottplib.LoadConfig(ussConfigFile,{'tolower':True})
if (ussCfg == None):
	ErrorExit("Error loading " + ussConfigFile)
	
logFile = os.path.join(logDir,'utcsteer.log') # this log will be common to several scripts
	
portSpeed = 9600
if ('hrog10:baud rate' in cfg):
	newSpeed = cfg['hrog10:baud rate']
	try:
		portSpeed = int(newSpeed)
	except:
		ErrorExit('Invalid baud rate:' + newSpeed)
	if (not portSpeed in BaudRates):
		ErrorExit('Invalid baud rate')

port = '/dev/hrog10'
if ('hrog10:port' in cfg):
	port = cfg['hrog10:port']

# Create the process lock		
lockFile = ottplib.MakeAbsoluteFilePath(cfg['main:lock file'],home,logDir)
Debug('Creating lock ' + lockFile)
if (not ottplib.CreateProcessLock(lockFile)):
	ErrorExit("Couldn't create a lock")

# Create UUCP lock for the serial port
uucpLockPath = '/var/lock';
if ('main:uucp lock' in cfg):
	uucpLockPath = cfg['main:uucp lock']
Debug('Creating uucp lock in ' + uucpLockPath)
ret = subprocess.check_output(['/usr/local/bin/lockport','-d',uucpLockPath,'-p',str(os.getpid()),port,sys.argv[0]])

if (re.match(rb'1',ret)==None):
	ottplib.RemoveProcessLock(lockFile)
	ErrorExit('Could not obtain a lock on ' + port + '.Exiting.')

Debug('Opening ' + port)

serport=None # so that this can flag failure to open the port
try:
	serport = serial.Serial(port,portSpeed,timeout=0.2)
except:
	Cleanup()
	ErrorExit('Failed to open ' + port)

serport.flush()

if args.show:
	[currTOFFS,currFFOF,currPHAS] = GetHROGSettings(serport)
	print('Time offset = {} ns'.format(currTOFFS))
	print('Fractional frequency offset = {}'.format(currFFOF))
	print('Phase offset = {} deg'.format(currPHAS))
	Cleanup()
	sys.exit(0)

UTCstr = 'UTC(' + UTCID +')'
html =  '<html>'
html += '<head></head>'
html += '<body>'
html += '<H2>' + UTCstr + ' steering advisory for '+  dt.strftime('%Y-%m-%d') + '</H2>'
html += '<br>'

html += 'Current MJD is ' + str(mjdToday) + '<br>'
html += '<br>'

# Check for a scheduled steer
steers = glob.glob(os.path.join(scheduledDir,'steer*.dat'))

newFFOF = None

if (len(steers)==0):
	Log(logFile,'no steers are scheduled')
	html += '<div> <strong> NO STEERS ARE SCHEDULED </div> </strong>'
elif (len(steers)==1):
	Debug('Found ' + steers[0])
	
	# Get current settings
	[currTOFFS,currFFOF,currPHAS] = GetHROGSettings(serport)
	msg = 'Current settings: TOFFS = {} FFOF = {} PHAS = {}'.format(currTOFFS,currFFOF,currPHAS)
	Log(logFile,msg)
	html += '<div> ' + msg + ' </div>'
	
	fin = open(steers[0],'r')
	
	for l in fin:
		if re.match(r'^\s*#',l):
			continue
		newFFOF = float(l)
		currFFOF = float(currFFOF)
		# Apply sanity checks to the frequency step
		# Never step more than a maximum amount and clamp the maximum frequency offset
		if (abs(newFFOF - currFFOF) > maxFSTEP) :
			msg = 'Scheduled steer {:e} exceeds permitted maximum step of {:e} - NOT STEERED'.format(newFFOF,maxFSTEP)
			Log(logFile,msg)
			html += '<div> <strong> ' + msg + ' </div> </strong>'
			newFFOF = None
		elif (abs(newFFOF) > maxFFOF):
			msg = 'Scheduled steer {:e} exceeds permitted maximum FFE of {:e} - NOT STEERED'.format(newFFOF,maxFFOF)
			Log(logFile,msg)
			html += '<div> <strong> ' + msg + ' </div> </strong>'
			newFFOF = None
			
		break
	fin.close()
else:
	Log(logFile,'too many unprocessed steers in {}'.format(scheduledDir))
	html += '<div> <strong> TOO MANY UNPROCESSED STEERS </div> </strong>'

if newFFOF:
	
	# FIXME check that values are as expected ?
	
	Debug('Steering')
	
	# Do the steer - this is total offset frequency ie the FFOF command
	if not(args.dryrun):
		HROGCmd(serport,'FFOF {:e}'.format(newFFOF)) # this gives 6 digits after the decimal point by default
	
	# note that the HROG will not necessarily echo back exactly what you sent eg 5.10 -> 5.1
	# Get new settings
	[toffs,ffof,phas] = GetHROGSettings(serport)
	msg = 'New settings: TOFFS = {} FFOF = {} PHAS = {}'.format(toffs,ffof,phas)
	Log(logFile,msg)
	html += '<div> ' + msg + ' </div>'.format(toffs,ffof,phas)
	
	# FIXME Check that the steer was nominal
	
	# Move the steering file
	if not(args.dryrun):
		shutil.move(steers[0],processedDir)
	
if email:
	
	msg = MIMEMultipart('related')
	msg['Subject'] = 'UTC steer  for ' + dt.strftime('%Y-%m-%d')
	msg['From'] = emailSender
	msg['To'] = recipients
	msg['Reply-To'] = emailSender

	body = MIMEText(html,'html')
	msg.attach(body)

	# Send the message via local SMTP server.
	s = smtplib.SMTP(SMTPserver)
	s.sendmail(emailSender,recipients.split(','), msg.as_string())
	s.quit()
	
Cleanup()
