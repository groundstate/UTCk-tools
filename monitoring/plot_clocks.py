#!/usr/bin/python3
# 
# Report on clocks contributing to UTC(AUS)
# The script will atempt to download UTCr data
# This allows plotting of UTC(AUS) and computing clock statistics wrt UTCr
# If the download fails, this is not a critical failure since UTC reporting must occur within a limited time window
# In this case, the plot of UTCr - UTC(AUS) is skipped and statistics are computed wrt UTC(AUS)

import allantools
import argparse
import datetime
import dateutil.relativedelta
import glob
import math

import numpy as np
import os
import re
import requests
import socket
import sys
import time

# This is where ottplib is installed
sys.path.append("/usr/local/lib/python3.6/site-packages")  # Ubuntu 18.04
sys.path.append("/usr/local/lib/python3.8/site-packages")  # Ubuntu 20.04
sys.path.append("/usr/local/lib/python3.10/site-packages") # Ubuntu 22.04
sys.path.append("/usr/local/lib/python3.12/site-packages") # Ubuntu 24.04

import ottplib as ottp

VERSION = "1.0.1"
AUTHORS = "Michael Wouters"

UTCR_LATENCY = 3     # in days
UTCR_REF_MJD = 56473 # Reference MJD for calculating end of week - last MJD for UTCr_1326

UTC_BASE_DIR = '/home/utc/utc'
CLK_DEFS_DIR = '/home/utc/utc/perl/process_clocks'

UTCID = 'AUS'

# ------------------------------------------
def Warn(msg):
	if (not args.nowarn):
		sys.stderr.write('WARNING! '+ msg+'\n')
	return

# ------------------------------------------

appName = os.path.basename(sys.argv[0])
debug = False

home =os.environ['HOME'] + '/'
user =os.environ['USER'] # remember to define this in the user's crontab

plotDir = './'
tmpDir  = './' 

utcBaseDir = UTC_BASE_DIR
clkDefinitionsDir = CLK_DEFS_DIR
historyLength = 3 # in months

bipmurl = 'https://webtai.bipm.org/api/v1.0'

parser = argparse.ArgumentParser(description='Plot clock data submitted to BIPM. The default is to report the previous month.',
	formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--debug','-d',help='debug (to stderr)',action='store_true')
parser.add_argument('--display',help='display plots',action='store_true')
parser.add_argument('--force',help='force UTCr download',action='store_true')
parser.add_argument('--history',help='history length (in months)',default = historyLength)
parser.add_argument('--month',help='month of report (1..12)')
parser.add_argument('--year',help='year  of report')
parser.add_argument('--utcbase',help = f'base directory for UTC reporting (default = {utcBaseDir})',default=utcBaseDir)
parser.add_argument('--tmp',help = f'directory for temporary files (default = {tmpDir} )',default = tmpDir)
parser.add_argument('--utc',help = 'generate file for UTC reporting',action='store_true')
parser.add_argument('--version','-v',action='version',version = appName + ' ' + VERSION + '\n' + 'Written by ' + AUTHORS)

args = parser.parse_args()

ottp.SetDebugging(args.debug)

if args.display:
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib.pyplot as plt
	import matplotlib.ticker as mticker
else:
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib as mplt # this (and the next line) stops warnings about being unable to connect to a display
	mplt.use('Agg') 
	import matplotlib.pyplot as plt
	import matplotlib.ticker as mticker

utcBaseDir = args.utcbase
tmpDir = args.tmp
historyLength = int(args.history)

if args.month and args.year:
	dt  = datetime.date( int(args.year),int(args.month),1) + dateutil.relativedelta.relativedelta(months=1) # pretend it's the next month
	tt = time.mktime(dt.timetuple())
	mjdToday = int(tt/86400)+40587
else:
	tt = time.time()
	mjdToday = int(tt/86400)+40587
	ottp.Debug(f'MJD today is {mjdToday:d}')

	dt = datetime.date.today() 
	dayToday = dt.day
	dt = dt.replace(day=1)

prevMonth = dt + dateutil.relativedelta.relativedelta(months=-1)
dstr = prevMonth.strftime('%Y_%b').lower()
fname = f'{dstr}.plots.all.pdf'
	
if args.utc:
	tmpDir  = os.path.join(home,'tmp') # keep it clean !
	plotDir = os.path.join(UTC_BASE_DIR,dstr) # it's a keeper
else:
	pass
	
prevMonth = dt + dateutil.relativedelta.relativedelta(months=-1)
dstr = prevMonth.strftime('%Y_%b').lower()
plotFileName = os.path.join(plotDir,f'{dstr}.plots.all.pdf')

# Get the names of clocks that we are reporting
# This is defined in clocks.* files

files = glob.glob(os.path.join(clkDefinitionsDir,'clocks.*'))
clocks = []
# clocks.append('UTC(Aus)') # this be a fudge
for f in files:
	ottp.Debug('Opening ' + f)
	fin = open(f,'r')
	headerRead = False
	while True:
		l = fin.readline()
		if not l: # EOF
			break
			
		l = l.strip()
		if not l: # skip blank lines
			continue
		
		if re.match('#ID1',l):
			headerRead = True
			continue
		
		if l[0] == '#': # skip comments
			continue
			
		if headerRead:
			ottp.Debug(f'Clock: {l}')
			clkLine = l.split()
			if len(clkLine) == 6:
				clocks.append(clkLine[0])
			else:
				ottp.Debug('Bad clock line')

	fin.close()

if args.debug:
	ottp.Debug('Clocks to be reported:')
	for c in clocks:
		ottp.Debug(c)
		
# UTCr is published each Wednesday 
# There's a 3 day lag (using rapid orbits and clocks for CGGTTS tweaking ?)
# Calculate the last expected MJD in UTCr
gotUTCr = False
nWeeks = int(math.floor((mjdToday - UTCR_REF_MJD)/7))
lastMJD = UTCR_REF_MJD + 7 * nWeeks
ottp.Debug('Last MJD in UTCr should be ' + str(lastMJD))

utcr=[[],[]]
# If the file is there, don't download (this is mainly for debugging)
futcr = os.path.join(tmpDir,'utcr.{:d}.dat'.format(lastMJD))
if (os.path.exists(futcr) and not(args.force)):
	ottp.Debug(f'UTCr file {futcr } exists - loading')
	fin = open(futcr,'r')
	for l in fin:
		if re.match('\s*#',l): # ignore comments
			continue
		ldata = l.split()
		if len(ldata) == 3: # v0.2 had two fields, v1.0 has three 
			try:
				utcr[0].append(int(ldata[0]))
				utcr[1].append(float(ldata[1]))
			except:
				ottp.Debug(f'Bad data in {fname}:{l}')
	fin.close()
	gotUTCr = True
else:
	ottp.Debug('Downloading UTCr')
		
	mjd1 = lastMJD - historyLength*31
	mjd2 = lastMJD

	httpreq = '{}/get-data.html?scale=utcr&lab={}&outfile=txt&mjd1={:d}&mjd2={:d}'.format(bipmurl,UTCID,mjd1,mjd2)
	try:
		resp = requests.get(httpreq)
		gotUTCr = True
	except:
		Warn('UTCr download failed')
		gotUTCr = False
	
	# Parse what we got back
	if gotUTCr:
		lines = resp.text.split('\r\n')
		for l in lines:
			if re.match('#',l): # ignore comments
				continue
			ldata = l.split()
			if len(ldata) == 3: # v0.2 of API had two fields, v1.0 has three 
				try:
					utcr[0].append(int(ldata[0]))
					utcr[1].append(float(ldata[1]))
				except:
					ottp.Debug('UTCRr bad data: ' + l)

		# Save the data for post-mortems
		fout = open(futcr,'w')
		fout.write(resp.text)
		fout.close()
		ottp.Debug('Saved UTCr data as ' + futcr)

if gotUTCr:
	firstUTCr = utcr[0][0]
	lastUTCr  = utcr[0][-1]
	ottp.Debug(f'UTCr data: first = {firstUTCr:d}, last = {lastUTCr:d}')

	# TOTDEV
	n = int(len(utcr[1])/2)
	taus = np.linspace(1,n,n+1)*86400
	(taus, totdevs, errors, ns) = allantools.totdev(utcr[1],rate = 1.0/86400.0,taus=taus)

	# TDEV
	n = int(len(utcr[1])/2)
	tdtaus = np.linspace(1,n,n+1)*86400
	(tdtaus, tddevs, errors, ns) = allantools.ttotdev(utcr[1],rate = 1.0/86400.0,taus=tdtaus)

	# FFE
	# Simple difference,  at midpoint
	# Don't use np.diff(), in case there is missing data
	fmjd = np.empty(len(utcr[0])-1)
	ffe  = np.empty(len(utcr[0])-1)

# Process clock data now that we have UTCr
clkData = {} # this will be used to store clock data 
for clk in clocks:
	clkData[clk]=[[],[],[]] # MJD , UTC(AUS) - clk, UTCr - CLK

for m in range(-historyLength,0):
	dtNew = dt + dateutil.relativedelta.relativedelta(months=m)
	dstr = dtNew.strftime('%Y_%b').lower()
	fname = os.path.join(utcBaseDir,dstr,f'{dstr}.all.diff')
	if os.path.exists(fname):
		fin = open(fname,'r')
		ottp.Debug(f'Opened {fname}')
	else: # try the archive directory
		fname = os.path.join(utcBaseDir,'YYYY_utc_archive',dtNew.strftime('%Y_utc'),dstr,f'{dstr}.all.diff')
		if os.path.exists(fname):
			fin = open(fname,'r')
			ottp.Debug(f'Opened {fname}')
		else:
			ottp.ErrorExit(f'Unable to open {fname}')
	
	# The files are strictly formatted but we won't assume too much
	headerRead = False
	clkCols = {}
	while True:
		l = fin.readline()
		if not l: # EOF
			break
			
		l = l.strip()
		if not l: # skip blank lines
			continue
		
		if re.match('Date',l):
			headerRead = True
			# The header contains clock identifiers
			# These can change month to month so we need to parse this line to identify the column to find clock data in
			header = l.split()
			ncols = len(header)
			# Clock names start in column 3
			for col in range(2,len(header)):
				clkCols[header[col]] = col
			continue
		
		if re.match('#page',l): # end of page 
			headerRead = False    # but there might be more
			continue
			
		if l[0] == '#': # skip comments TODO CTRL-L in multi-page files
			continue
			
		# Should have something parseable now
		data = l.split()
		if len(data) == ncols:
			mjd = int(data[1])
			if gotUTCr:
				utcrdiff = None
				if mjd >= firstUTCr and mjd <= lastUTCr:
					# Check the MJD
					imjd = mjd - firstUTCr
					if utcr[0][imjd] == mjd:
						utcrdiff = utcr[1][imjd]
					else:
						Warn(f'Failed to find MJD {mjd:d} in UTCr file') # shouldn't happen 
			else:
				utcrdiff = 0.0 # this is so we can compute clkdiff + utcrdiff
			for clk in clocks:
				if clk in clkCols:
					col = clkCols[clk]
					if ('*' in data[col]):
						ottp.Debug(f'Missing data line in {l}')
					else:
						clkdiff = float(data[col])
						clkData[clk][0].append(mjd)
						clkData[clk][1].append(clkdiff)
						if None != utcrdiff:
							clkData[clk][2].append(clkdiff + utcrdiff) # 
						else:
							clkData[clk][2].append(None)
		else:
			ottp.Debug(f'Bad data line {l}')
	
	fin.close()

# The list clkData will almost certainly have missing UTCr data at the end. 
# These are flagged by clkData[2][] == None
# Make a new list that has missing points removed so that we can do statistics on valid data

cleanClkData = {} # this will be used to store clock data 
for clk in clocks:
	cleanClkData[clk]=[[],[],[]]
	
for clk in clocks:
	i = 0
	while i < len(clkData[clk][0]):
		if (clkData[clk][2][i] != None):
			cleanClkData[clk][0].append(clkData[clk][0][i])
			cleanClkData[clk][1].append(clkData[clk][1][i])
			cleanClkData[clk][2].append(clkData[clk][2][i])
		i = i + 1
			
if gotUTCr:
	for i in range(1,len(utcr[0])):
		fmjd[i-1] =(utcr[0][i] + utcr[0][i-1])/2.0
		ffe[i-1] = -1.0E-9*(utcr[1][i] - utcr[1][i-1]) / ((utcr[0][i] - utcr[0][i-1])*86400) # correct units Hz/Hz


# Now produce some plots
	
with PdfPages(plotFileName) as pdf:
	
	# SECTION 1
	# Plot UTCr-UTC(AUS)
	if gotUTCr:
		fig,(ax1,ax2,ax3,ax4) = plt.subplots(4,1,figsize=[8,12])
		dstr = prevMonth.strftime('%B %Y')
		title = 'Clock performance for ' + dstr + '\n'
		title += plotFileName + '\n'
		title += os.path.basename(sys.argv[0])+ ' v' + VERSION   + '     run ' + datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S') + '\n'
		fig.suptitle(title,ha='left',x=0.1)
				
		ax1.plot(utcr[0],utcr[1])
		ax1.scatter(utcr[0],utcr[1],s=9,marker='o')
		ax1.set_ylabel('delta (ns)')
		ax1.set_xlabel('MJD')
		ax1.set_title('UTCr - UTC(AUS)')
		ax1.grid()

		ax2.plot(fmjd,ffe/1.0E-14)
		ax2.scatter(fmjd,ffe/1.0E-14,s=9,marker='o')
		ax2.set_ylabel(r'ffe ($10^{-14}$ Hz/Hz)')
		ax2.set_xlabel('MJD')
		ax2.grid()

		ax3.loglog(taus/86400,totdevs*1.0E-9/1.0E-14,'o-')
		ax3.xaxis.set_major_formatter(mticker.ScalarFormatter())
		# ax3.xaxis.set_minor_formatter(mticker.ScalarFormatter())
		ax3.yaxis.set_major_formatter(mticker.ScalarFormatter())
		ax3.yaxis.set_minor_formatter(mticker.ScalarFormatter())
		ax3.set_ylabel(r'frac. TOT DEV ($10^{-14}$ Hz/Hz)')
		ax3.set_xlabel('tau (days)')
		ax3.grid(visible=True,which='both',axis='both')

		ax4.loglog(tdtaus/86400,tddevs,'o-')
		ax4.xaxis.set_major_formatter(mticker.ScalarFormatter())
		ax4.yaxis.set_minor_formatter(mticker.ScalarFormatter())
		ax4.set_ylabel('TIME TOT DEV (ns)')
		ax4.set_xlabel('tau (days)')
		ax4.grid(visible=True,which='both',axis='both')
		
		pdf.savefig()
	
	# SECTION 2
	# We plot the clocks using the data that are reported to BIPM
	# since we want a visual check on what has been assembled for submission
			
	newPage = True
	pageComplete = False
	iPlot = 0
	for clk in clocks:
		
		if not clkData[clk][0]:
			continue
		
		iPlot += 1
		if iPlot == 7:
			newPage = True
			iPlot = 1
			pageComplete = True
			
		if newPage:

			if pageComplete:
				pageComplete = False
				pdf.savefig()
				
			fig,axs = plt.subplots(3,2,figsize=[8,12])	
			
			dstr = prevMonth.strftime('%B %Y')
			if not gotUTCr:
				title = 'Clock performance for ' + dstr + '\n'
				title += fname + '\n'
				title += 'WARNING! Failed to download UTCr. All statistics computed wrt UTC(AUS)\n\n'
				title += os.path.basename(sys.argv[0])+ ' v' + VERSION   + '     run ' + datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S') + '\n'
				fig.suptitle(title,ha='left',x=0.1)
			else:
				fig.suptitle('Clock history')
			newPage = False
		prow = int(math.floor((iPlot-1) / 2))
		pcol = (iPlot-1) - prow*2
		axs[prow,pcol].plot(clkData[clk][0],clkData[clk][1],label='UTC(AUS) - clk')
		if gotUTCr:
			axs[prow,pcol].plot(clkData[clk][0],clkData[clk][2],label='UTCr - clk')
		#plt.scatter(clkData[clk][0],clkData[clk][1],s=9,marker='o')
		axs[prow,pcol].grid()
		axs[prow,pcol].legend(fontsize=6)
		axs[prow,pcol].set_title(clk)
		
	pdf.savefig() # catch the last one
	
	# SECTION 3
	# Time history again, but with rate and offset subtracted
	newPage = True
	pageComplete = False
	iPlot = 0
	for clk in clocks:
		
		if not cleanClkData[clk][0]:
			continue
		
		iPlot += 1
		if iPlot == 7:
			newPage = True
			iPlot = 1
			pageComplete = True
			
		if newPage:

			if pageComplete:
				pageComplete = False
				pdf.savefig()
				
			fig,axs = plt.subplots(3,2,figsize=[8,12])	
			fig.suptitle('Clock history - rate and offset removed')
			newPage = False
		
		tclk = np.array(cleanClkData[clk][0])
		uclk = np.array(cleanClkData[clk][1])
		coeff = np.polyfit(tclk,uclk,1) # units are ns/day
		foffset = -coeff[0] # note the sign!
		toffset =  coeff[1]
		
		prow = int(math.floor((iPlot-1) / 2))
		pcol = (iPlot-1) - prow*2
		axs[prow,pcol].plot(tclk,uclk + foffset*tclk - toffset,label='UTC(AUS) - clk')
		if gotUTCr:
			uclk = np.array(cleanClkData[clk][2])
			coeff = np.polyfit(tclk,uclk,1) 
			foffset = -coeff[0] # note the sign!
			toffset =  coeff[1]
			axs[prow,pcol].plot(tclk,uclk + foffset*tclk - toffset,label='UTCr - clk')
		axs[prow,pcol].grid()
		axs[prow,pcol].legend(fontsize=6)
		axs[prow,pcol].set_title(clk)
		
	pdf.savefig() # catch the last one

	# SECTION 4
	# Clock stability
	fig,(ax1,ax2) = plt.subplots(2,1,figsize=[8,12])	
	fig.suptitle('Clock stability')
	
	for clk in clocks:
		if not clkData[clk][0]:
			continue
		dclean = [d for d in clkData[clk][2] if d!= None] # remove missing points
		n = int(len(dclean)/2)
		taus = np.linspace(1,n,n+1)*86400
		(taus, totdevs, errors, ns) = allantools.totdev(dclean,rate = 1.0/86400.0,taus=taus)
		ax1.loglog(taus/86400,totdevs*1.0E-9/1.0E-14,label=clk)
	
	ax1.legend(fontsize=6)
	ax1.xaxis.set_major_formatter(mticker.ScalarFormatter())
	ax1.yaxis.set_minor_formatter(mticker.ScalarFormatter())
	ax1.yaxis.set_major_formatter(mticker.ScalarFormatter())
	ax1.set_ylabel(r'frac. TOT DEV ($10^{-14}$ Hz/Hz)')
	ax1.set_xlabel('tau (days)')
	ax1.grid(visible=True,which='both',axis='both')
	
	for clk in clocks:
		if not clkData[clk][0]:
			continue
		dclean = [d for d in clkData[clk][2] if d!= None] # remove missing points
		n = int(len(dclean)/2)
		tdtaus = np.linspace(1,n,n+1)*86400
		(tdtaus, tddevs, errors, ns) = allantools.ttotdev(dclean,rate = 1.0/86400.0,taus=tdtaus)
		ax2.loglog(tdtaus/86400,tddevs,label=clk)
	
	ax2.legend(fontsize=6)
	ax2.xaxis.set_major_formatter(mticker.ScalarFormatter())
	ax2.yaxis.set_minor_formatter(mticker.ScalarFormatter())
	ax2.yaxis.set_major_formatter(mticker.ScalarFormatter())
	ax2.set_ylabel('TIME TOT DEV (ns)')
	ax2.set_xlabel('tau (days)')
	ax2.grid()
	ax2.grid(visible=True,which='both',axis='both')
	
	pdf.savefig()
	
# Display all plots when we're done
if args.display:
	plt.show()

# Clean up any temporary files
if not args.debug:
	os.unlink(futcr)
print(f'Plots saved to {plotFileName}')
