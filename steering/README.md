# Description
Steering is based on feedback from UTCr, which is updated each week, late on Wednesdays (UTC).

The steering is done is two steps.

First, the script utcsteersched.py is run periodically.
It calculates the current MJD and decides whether new UTCr data might be available to download.
If it's too early for new data to be available, the script exits.
If it's more than 3 days after UTCr should have been published, the script exits.
If the data have already been downloaded (determined by checking a status file), the script exits.
These three checks are bypassed if the --force option was used.

Download of new data is then attempted.
If this is successful then the steering parameters are computed and stored in a file.
Information about the steer is then collated into an email and sent out for checking.

The second step is to implement the steer.
The script  utcsteer.py does this. 
It is written for a Spectradynamics HROG.
The script is meant to be run a few hours after utcsteersched.py, with the window coinciding within
working hours. This is to allow time to examine the steer, and intervene if necessary.

# Installation

You need

- numpy  (python3-numpy)  
- scipy  (python3-scipy)
- matplotlib (python3-matplotlib)
- htmldoc (htmldoc,htmldoc-common)
- allantools (pip3 install allantools)

Run install.py to create the necessary directories etc.


# Using with cron

You may need to define USER and proxy servers in your crontab
eg in the crontab

USER=utcsteer


# Configuration file

The sample configuration file documents the available configuration options.

