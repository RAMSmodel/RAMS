#!/bin/bash
################################################################################
# This script will run a single test serially (by default) (on any system) or in
# parallel (on a simple computer cluster that does not use a PBS/QSUB queuing
# system). Note that many supercomputing systems (that use modules and PBS
# queuing systems) often have different methods for compiling code for serial
# or parallel processing. The correct compiling method is required for running
# this script in parallel or serial.
#
# A couple of important items here for this test simulation.
# 1. RAMSIN.testrunonly - produces a 3D 1-grid supercell simulation with
#  LEVEL=3 microphysics
# 2. This idealized simulation does not require geographical data (sfctypehdf5)
#  for running since geography does not matter. Also, we do not need gridded
#  initialization data beyond the sounding input in the namelist. The sounding
#  is applied horizontally homogeneous for an idealized run.
#
# Here is a sample executable statement for running RAMS in parallel,
#  but NOT on a supercomputer using a PBS QSUB queuing system. Consult their
#  userguides for running parallel jobs on their systems. Each system is unique.
# /home/smsaleeb/software/mpich-3.3.2/bin/mpiexec -machinefile machs -np 8 \
# ./bin.rams/rams-6.3.01 -f RAMSIN.testrunonly
# Note: might have to add (-iface eth0) to this executable statement or
#  something similar depending on how your compute nodes communicate. Use
#  utility (ifconfig) to find out how your nodes communicate (ie. eth0 or eth1).
#
# Do not blindly use this script for running every simulation. It is simply
# included here as a starting point and because I use this frequently on my
# desktop Linux workstation for RAMS tests and debugging. If using this script
# and RAMSIN namelist, please become familiar with its flags, settings, and
# functionality. They cannot be used for every runtime situation. Many RAMSIN
# flags in this test example are turned on for testing only and are not
# appropriate for a real simulation to study.
###############################################################################
# Set your RAMS root path
rd=`pwd`/..
# RAMSIN name
ramsin="RAMSIN.testrunonly"
# RAMS version (ie. 6.1.6)
vs=6.3.04
# Set flag for type of test (0=sequential, 1=parallel)
runtype=1
# Set number of nodes for parallel run.
n=8
# Set delete flag (0 = do not delete, 1 = delete and start over)
del=1

# Flag for whether to create a machines file "machs" or not. (0=no, 1=yes)
# Set machsmake=0 if you have a custom file. Set machsmake=1 to create one
# based on Hostname of this local machine.
machsmake=1

# Set your parallel executable here
a1=/home/smsaleeb/software/mpich-3.3.2/bin/mpiexec
# Set you machines file path here
a2=$rd/bin.rams/machs
# Set RAMS executable path and name
a3=$rd/bin.rams/rams-$vs
# Set REVU post-processor executable path and name
a4=$rd/bin.revu/revu-$vs

#Check to see that RAMSIN namelist exists
if [ ! -f $rd/bin.rams/$ramsin ]; then
 echo "Input paths for this file is incorrect: RAMSIN.testrunonly"
 exit
fi

#Check to see that RAMS executable exists
if [ ! -f $a3 ]; then
 echo "Cannot find RAMS executable such as: bin.rams/rams-6.3.02"
 echo "Need to compile RAMS for this script to work."
 exit
fi

#If trying a parallel run, make sure machines file exists.
if [ $runtype -eq 1 ]; then
  if [ $machsmake -eq 0 ]; then
    if [ ! -f $a2 ]; then echo "You need a machines file 'machs'"; exit; fi;
  elif [ $machsmake -eq 1 ]; then
    echo $HOSTNAME > $a2
  else
    echo "Need to set MACHSMAKE = 0 or 1"
    exit
  fi
fi

#Check to see that MPI executable and machines files exist for a parallel test
if [ $runtype -eq 1 ]; then
 if [ ! -f $a1 ] || [ ! -f $a2 ]; then
  echo "One of your parallel simulation input paths is incorrect. Stopping!"
  echo "MPI executable or machs file"
  exit
 fi
fi

# DONE WITH NECESSARY USER CHANGES
###############################################################################

for dirname in $rd/bin.rams/testrun.output
do
  if [ -d $dirname -a $del -eq 1 ]; then
   rm -f $dirname/*.h5 $dirname/*.txt
  else
   if [ ! -d $dirname ]; then mkdir $dirname; fi;
  fi
done

rc="$a3 -f"

#Check on runtype settings for sequential or parallel
if [ $runtype -eq 0 ]; then
  rc1=$rc
elif [ $runtype -eq 1 ]; then
  rc1="$a1 -machinefile $a2 -np $n $a3 -f"
else
  echo "Need to set RUNTYPE = 0 or 1"
  exit
fi

$rc1 $rd/bin.rams/$ramsin
