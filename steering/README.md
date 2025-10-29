# Installation

You need

- numpy  (python3-numpy)  
- scipy  (python3-scipy)
- matplotlib (python3-matplotlib)
- htmldoc (htmldoc)
- htmldoc  (htmldoc-common)
- allantools (pip3 install allantools)

Run install.py to create the necessary directories etc


# Using with cron

You will likely need to define USER and proxy servers in your crontab
eg in the crontab

USER=utcsteer


# Configuration file

## [Main] section 
### history

The number of days that are fetched for UTCr

### email recipients

This is a comma-separated list

### reports

directory for HTML reports

### log 

directory for log files

### control 

directory for control and status files

Within this directory are the directories processed_steer and scheduled_steer
Pending steer files are written to scheduled_steer
When they are successfully processed, they are moved to scheduled_steer

### tmp

directory for temporary files like UTCr data
