#!/bin/sh

# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=$(dirname "$SCRIPT")
PROJECTNAME=$(basename "$SCRIPTPATH")

BTYPE=RelWithDebInfo
BIN=${SCRIPTPATH}-build/${BTYPE}/${PROJECTNAME}

${BIN} plot_NBrS_XS_exact_Ar.xml | tee Log.txt 
