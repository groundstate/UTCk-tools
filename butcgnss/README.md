# Description

This is a set of scripts that can be used to generate reports on the differences UTC(k)-bUTC_GNSS and UTC-bUTC_GNSS.
In particular, they are intended to provide a data flow suitable to the requirements of the quality system that an NMI
will typically be operating under. The documentation here assumes that you are installing on Linux.

The scripts define the following data flow:
- butcupdate.py  Updates the database of UTC(k)-bUTC_GNSS and UTC-bUTC_GNSS differences (automatic)
- butcreport.py  Updates the monthly reports from the database (automatic)
- butccheck.py   Produce a report on computed differences for examination (automatic)
- butcrelease.py Release data for final publication (manual)

# Installation

Most of the python modules you need will likely already be installed.

You will probably need to install cggttslib, ottplib and rinexlib though.
These are available through OpenTTP:
https://github.com/openttp/openttp

Check out 'develop', cd to software/system/ and then, with superuser privileges,
./installys.py -i cggttslib
./installys.py -i ottplib
./installys.py -i rinexlib

# Configuration file

A single configuration file controls the operation of all four scripts.
The sample configuration file documents the available configuration options.


# Setting up cron

The first three scripts listed above are intended to be run automatically.

Since the differences UTC(k) - bUTC_GNSS can be computed each day, butcupdate.py can be run daily.
The caveat on this is that uB is needed for the calculation of uncertainties and this is not strictly 
available until the publication of Circular T.

The current software uses the last value of uB published in Circular T, when the current value is not obtainable. 
In practice, the incremental change in the published uncertainty due to ageing of the calibration will at
most affect the published uncertainty by +/- 1 ns.
 
butcreport.py can be run each day. It recalculates all values in the database within a window of XX days, 
from the current day.

butcrcheck.py is run after Circular T is published.

butcrelease.py is run manually.
It indicates the MJDs for which data updated from Circular T is now available.
Values of UTC-bUTC_GNSS are not reported in he published report until they are released.


You may need to define USER in your crontab
eg in the crontab
USER=butcgnss






