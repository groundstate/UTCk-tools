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
import json
import os
import requests
import sys

try:
	import ottplib   as ottp
except ImportError:
	sys.exit('ERROR: Must install ottplib\n eg openttp/software/system/installsys.py -i ottplib')

VERSION = '0.0.1'
AUTHORS = 'Michael Wouters'

# ------------------------------------------
def ShowVersion():
	print (os.path.basename(sys.argv[0])+" "+VERSION)
	print ('Written by ' + AUTHORS)
	return

# ---------------------------------------------
def FetchBUTC(gnss,startMJD,stopMJD):
	ottp.Debug(f'Fetching Circular bUTC prediction data  for {gnss} from BIPM');
	req = f'{httpRequest}scale=b_{gnss}&mjd1={startMJD}&mjd2={stopMJD}&outfile=json'
	ottp.Debug(req)
	try:
		r = requests.get(req,verify=rootCert)
	except:
		ottp.Debug('http request failed')
		return None,None,None
	
	d = json.loads(r.text)
	if not(d['errorcode'] == 0):
		ottp.Debug(f'http request return error {d['errorcode']}')
		return None,None,None
	
	#data = [[m,b,u] for m,b,u in zip(d['data'][0]['x'],d['data'][0]['y'],d['data'][0]['unc'])]

	return d['data'][0]['x'],d['data'][0]['y'],d['data'][0]['unc']
	
# ----------------------------------------------------------------------------------
home =os.environ['HOME']
root = home

httpRequest = 'https://webtai.bipm.org/api/v1.0/get-data.html?'
rootCert = None # if you use an empty string, this will skip SSL verification which is a bad thing

configFile = os.path.join(root,'etc/butc.conf')
tmpDir  = os.path.join(root,'tmp')

if ottp.LibMajorVersion() >= 0 and ottp.LibMinorVersion() < 2: 
	sys.exit('Need ottplib minor version >= 2')

examples='TO DO'
parser = argparse.ArgumentParser(description='Provides information for checking a monthly report',
	formatter_class=argparse.RawDescriptionHelpFormatter,epilog=examples)
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

cfg=ottp.Initialise(configFile,['main:gnss'])

gnss = cfg['main:gnss'].split(',')
gnss = [g.strip() for g in gnss] 

if ('main:root certificate' in cfg):
	rootCert= cfg['main:root certificate']

# Make a plot for each GNSS
for g in gnss:
	m,d,u = FetchBUTC(g,60100,60110)  # nb u is string!
	
	
