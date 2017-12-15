#!/bin/bash
# Linux script to support clean calling of Matlab as a 
# command processor for use in make scripts.
#
# Usage:
# >> matlabshell <MATLAB_SCRIPT>
#
# Installation: 
# - The first line must point to a valid shell binary. 
#   This script originally written for the Bash shell.
#
# - Matlab must be in the shell search path.
#
# - LOGFILE must point to a valid directory.
#
# Doug Harriman (http://www.linkedin.com/in/dougharriman)
#

LOGFILE=~/tmp/matlab-shell.log
rm -f $LOGFILE
stty echo
matlab -nodisplay -nosplash -r CC_clusterer_BNCfeaturesGW_noun.m
# Return with Matlab's return code
exit $?