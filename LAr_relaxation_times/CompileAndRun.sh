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
rm -rf ${SCRIPTPATH}/${BTYPE}
mkdir ${SCRIPTPATH}/${BTYPE}
cd ${SCRIPTPATH}/${BTYPE}
# set -x displays cmake command. Brackets create subshell so that there is no need ot call set +x
(set -x; cmake -DCMAKE_BUILD_TYPE=${BTYPE} -DBOOST_ROOT=${BOOSTPATH} ../)
make
make target install
rm -f ../Log.txt
./${PROJECTNAME} ../settings.xml | tee ../Log.txt

