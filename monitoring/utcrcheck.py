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
# This should run ...

import argparse
import hashlib
import math
import ftplib
import io
import os
from pathlib import Path
import sys
import time

# This is where ottplib is installed
sys.path.append("/usr/local/lib/python3.6/site-packages")  # Ubuntu 18.04
sys.path.append("/usr/local/lib/python3.8/site-packages")  # Ubuntu 20.04
sys.path.append("/usr/local/lib/python3.10/site-packages") # Ubuntu 22.04
sys.path.append("/usr/local/lib/python3.12/site-packages") # Ubuntu 24.04

import ottplib as ottp

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
			
VERSION = "0.1.0"
AUTHORS = "Michael Wouters"

DELAY  = 1  # delay between current MJD and most recent upload
NCHECK = 2  # number of MJDs to go back for check

home =os.environ['HOME'] + '/'
user =os.environ['USER'] # remember to define this in the crontab
configFile = os.path.join(home,'etc/utcrcheck.conf')
tmpDir = os.path.join(home,'tmp')

tt = time.time()
mjdToday = int(tt/86400)+40587

parser = argparse.ArgumentParser(description='Check clock data uploaded to BIPM',
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

# Useful for debugging problems
# ftp = ftplib.FTP(server)
# ftp.set_debuglevel(3)
# ftp.sanitize = lambda s: repr(s) # makes password visible
# ftp.login(user=usr, passwd=password)

# Construct the list of files to download and check
files = []
for m in range(mjdToday - DELAY- NCHECK + 1, mjdToday- DELAY + 1):
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
		for i in range(0,2):
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

if checkOK:
	for f in files:
		rep = f[0]
		if not(f[1]):
			rep += ' local file is missing'
		elif not(f[2]):
			rep += ' remote file is missing'
		elif not(f[3]):
			rep += ' bad checksum'
		else:
			rep += ' OK'
		print(rep)
else;
	pass
