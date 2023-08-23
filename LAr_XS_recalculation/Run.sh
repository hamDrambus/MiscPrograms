#!/bin/sh

# Necessary input for CMake
#GEANT4PATH=${HOME}/Software/Geant4/geant4-v11.0.0-install/lib/Geant4-11.0.0/
BOOSTPATH=${HOME}/Software/boost_1_67_0/

# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=$(dirname "$SCRIPT")
PROJECTNAME=$(basename "$SCRIPTPATH")

#-DCMAKE_BUILD_TYPE=RelWithDebInfo
#-DCMAKE_BUILD_TYPE=Debug
BTYPE=RelWithDebInfo
# Clear build directory and cd to it
cd ${SCRIPTPATH}/${BTYPE}
rm -f ../Log.txt
./${PROJECTNAME} ../settings.xml | tee ../Log.txt

