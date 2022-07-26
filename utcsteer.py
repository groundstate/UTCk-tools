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
import subprocess
import sys
import time


# This is where ottplib is installed
sys.path.append("/usr/local/lib/python3.6/site-packages") # Ubuntu 18.04
sys.path.append("/usr/local/lib/python3.8/site-packages") # Ubuntu 20.04

import ottplib

VERSION = "0.0.1"
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
def SendHROGCmd(serport,cmd):
	Debug('Sending ' + cmd)
	ba = bytearray() # stupid unicode
	cmd = cmd + '\r'
	ba.extend(cmd.encode())
	serport.write(ba)
	
#-----------------------------------------------------------------------------
def GetHROGResponse(serport,cmd):
	
	for i in range(0,5):
		SendHROGCmd(serport,cmd)
		d = serport.readline().decode('utf-8') # assuming  here that ECHO is on 
		if not(re.search(cmd,d)):
			continue
		d = serport.readline().decode('utf-8').strip()
		if re.search(cmd,d):
			Debug('Received ' + d)
			return d.split() # return a list of strings since usual response is CMD VALUE UNIT 
	return []

#-----------------------------------------------------------------------------
def GetHROGSettings(serport):
	stoffs = None
	ffof = None
	phase = None
	ret = GetHROGResponse(serport,'STOFFS?') # time offset STOFFS?
	if ret:
		stoffs = float(ret[1]) # WARNING assuming units are always ns
	ret = GetHROGResponse(serport,'FFOF?') # current fractional frequency offset FFOF?
	if ret:
		ffof = float(ret[1])
	ret = GetHROGResponse(serport,'PHAS?')# phase offset PHAS?
	if ret:
		 phas = float(ret[1])
	return [stoffs,ffof,phas]

# ------------------------------------------
# Main

appName = os.path.basename(sys.argv[0])
debug = False

home =os.environ['HOME'] + '/'
user =os.environ['USER'] # remember to define this in the user's crontab
configFile = os.path.join(home,'etc/utcsteer.conf')
repDir  = os.path.join(home,'utcsteer/reports')
logDir = os.path.join(home,'utcsteer/logs')
controlDir = os.path.join(home,'utcsteer/control')
scheduledDir = os.path.join(controlDir,'scheduled_steer')
processedDir = os.path.join(controlDir,'processed_steers')

parser = argparse.ArgumentParser(description='Executes a steer calculated and scheduled by utcsteersched.py',
	formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--config','-c',help='use an alternate configuration file',default=configFile)
parser.add_argument('--debug','-d',help='debug',action='store_true')
parser.add_argument('--show','-s',help='show HROG settings',action='store_true')
parser.add_argument('--version','-v',action='version',version = os.path.basename(sys.argv[0])+ ' ' + VERSION + '\n' + 'Written by ' + AUTHORS)

args = parser.parse_args()

debug = args.debug

configFile = args.config;

if (not os.path.isfile(configFile)):
	ErrorExit(configFile + ' not found')
	
logPath = os.path.join(home,'logs')
if (not os.path.isdir(logPath)):
	ErrorExit(logPath + "not found")

cfg=Initialise(configFile)

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
	[stoffs,ffof,phas] = GetHROGSettings(serport)
	print('Time offset = {:f} ns'.format(stoffs))
	print('Fractional frequency offset = {:f}'.format(ffof))
	print('Phase offset = {:f} deg'.format(phas))
	Cleanup()
	sys.exit(0)
