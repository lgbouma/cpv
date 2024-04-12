#!/bin/bash

# Description:
#   This script monitors the status of the job "python -u find_CPVs.py" and
#   sends an email notification to the email address specified in the file
#   "~/.secret_email" if the job is not running. The script checks the status
#   of the job every 60 seconds and exits after sending the first email.
#
# Usage:
#   ./monitor_CPV_finder.sh &
#
# Note:
#   - Make sure to have 'ssmtp' installed and properly configured on your
#     system to send emails from the command line.
#   - Create a file named "~/.secret_email" and add your email address
#     (e.g., whatever@gmail.com) as a single line in the file.

# Read the email address from the file
email_address=$(cat ~/.secret_email)

while true; do
	  if ! pgrep -f "python -u find_CPVs.py" > /dev/null; then
        echo -e "To: $email_address\nSubject: ALERT: find_CPV not running\n\n"\
"The job 'python -u find_CPVs.py' is not running." | ssmtp $email_address
        break
    fi
    sleep 60
done
