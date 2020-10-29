#! /bin/bash
##############################################################
# Program: x.repack.sh
# Programmer: Robert Seigel
#             University of Miami
#             rseigel@rsmas.miami.edu
#             16 March 2013
# Execute: x.repack.sh <type> [<ramssubprocess>] [<pattern>]
#          <type>: surface or analysis MANDATORY
#          [<ramssubprocess>]: MANDATORY if <type> = analysis
#                              If using LSF job execution, then 
#                              <ramsubprocess> = Job Name
#                              If NOT using LSF job execution, then 
#                              <ramsubprocess> = 6.0-opt, if the
#                              executable is rams-6.0-opt
#          [<pattern>]: OPTIONAL 3rd argument if RAMS has already
#                       finished and the user wants to repack all 
#                       files in the supplied directory
#                       *** dummy argument 2 still required ***
# Purpose: This script repacks the HDF5 files while a distributive
#          memory RAMS run is occurring. Because Parallel HDF5 is
#          not yet supporting online compression, this script 
#          takes its place to conserve disk space.
# Overview: This script runs while RAMS is running. It first looks
#           in RAMSIN to find where the files are being written.
#           Then, it iteratively checks to see if new files are
#           created. If a new file is output, it waits until the
#           file is created and then uses h5repack to repack
#           the HDF5 file. Once the RAMS processes are finished,
#           x.repack.sh will terminate.
# 
#******************************************************************
# It is suggested to run this process in the background, as 
# simulations often take a long time. To do so, here is an example:
#    nohup x.repack.sh analysis 6.0-opt > repack.txt &
# This will execute x.repack.sh to repack the HDF5 analysis files
# that are output from rams-6.0-opt. The script output is placed
# in repack.txt for reference.
##############################################################


##############################################################
# FUNCTION to find the new file and repack it according to 
# specified gzip level
repack ()
{
   echo "Start repack"

   # Initialize some variables if this is our first entry into this routine
   if [ $lastmod -eq 0 ]; then
      newmod=0
      oldfsize=0
   fi

   # Loop over .h5 files in requested directory
   for f in $thisdir*'.h5'
   do

      # get time since file was last modified
      modsecs=$(date --utc --reference=$f +%s)
      nowsecs=$(date +%s)
      delta=$(($nowsecs-$modsecs))
      filesize=$(stat -c%s "$f")

      # Test if the following are true in order to repack the file
      # (1) Current file AGE is LE user-defined age threshold (filter out old sims)
      # (2) Current file is EITHER (1) newer than previous file OR
      #                            (2) larger than previous file
      # (3) Current file SIZE is GT user-defined size threshold (filter out already
      #                                                          repacked files)  
      if [[ $delta -le $repacksecs ]] && \
	 [[  ( $modsecs -gt $lastmod || $filesize -gt $prevfilesize ) ]] && \
         [[  $filesize -gt $unpackFsize ]]; then

         # Wait UNTIL header file exists to make sure HDF5 write is done
         if [ $runtype == "analysis" ]; then
            hfile=${f:0:$((${#f}-5))}'head.txt'

            until [ -e $hfile ]; do sleep 1; done

            # Case where file is being overwritten
            until [ ! $hfile -ot $f ]; do sleep 2; done
   
            # In case the header file is big
            sleep 3
         fi

         # Repack the file
         echo "h5repack -f SHUF -f GZIP=$gziplevel $f $f.temp"   #***UNCOMMENT
         h5repack -f SHUF -f GZIP=$gziplevel $f $f.temp   #***UNCOMMENT
         mv $f.temp $f                                    #***UNCOMMENT
         echo "New file $f has been repacked because it was created $delta secs ago"
         newmod=$(date --utc --reference=$f +%s)

         # Now copy the file to a remote computer. RSA keys have already been generated.
	 #scp $f ldgrant@frost.atmos.colostate.edu:$supcelldir
	 #scp $hfile ldgrant@frost.atmos.colostate.edu:$supcelldir
#         scp $f $copydir 
#         scp $hfile $copydir 
                  
      fi
      prevfilesize=$(stat -c%s "$f")
   
   done
   lastmod=$newmod
}
######################################################################################
######################################################################################
# ONLYEDIT THESE

#set -x

# Define some parameters
exectype=1         # Type of job execution.
                   # 1 = LSF
                   # 2 = Standard way, using a direct call to MPI
runtype=$1         # this is either surface or analysis from command line
ramsname=$2        # second command line input *MANDATORY for analysis*
                   # This is not used if exectype=1
#repacksecs=108000  # maximum age in seconds of a file that is to be repacked (30 h)
repacksecs=1080000  # maximum age in seconds of a file that is to be repacked (30 h)
unpackFsize=5000  # megabyte threshold for unpacked file size (~ 5G)
#unpackFsize=1  # megabyte threshold for unpacked file size (~ 1M)
gziplevel=6        # gzip level
checkint=120       # in seconds
#supcelldir='/tasman/ldgrant/LPsup_300/' # Directory for file transfers
#copydir='ldgrant@frost.atmos.colostate.edu:/tasman/ldgrant/LPsup_300/z.'$2'/' # Directory for file transfers
#echo "copydir: $copydir"

# END ONLYEDIT
######################################################################################
######################################################################################
# Perform some checks

# Did the user input the process name during an analysis repack?
if [ "$ramsname" == "rams-" -a "$runtype" == "analysis" ]; then
   echo 'Please specify the RAMS executable name, excluding rams-'
   exit
fi
# Is the user exercising the 3rd argument? If so, dont read RAMSIN
if [ -z "$3" ]; then
   RAMSINread=1
else
   RAMSINread=0
fi
# Convert file size to megabytes
unpackFsize=$((unpackFsize*1024*1024))

##############################################################
# First Read RAMSIN to grab file locations
# This method is based on the fact that each RAMSIN parameter
# is preceeded by THREE SPACES. 
if [ $RAMSINread -eq 1 ]; then
   # TOPFILES
   # find the line
   topdir=`egrep '^   TOPFILES' RAMSIN`
   # remove suffix from 2nd single quote
   topdir=${topdir%\'*}
   # remove prefix from 1st single quote
   topdir=${topdir#*\'}
   
   # SFCFILES
   sfcdir=`egrep '^   SFCFILES' RAMSIN`
   sfcdir=${sfcdir%\'*}
   sfcdir=${sfcdir#*\'}
   
   # SSTFILES
   sstdir=`egrep '^   SSTFPFX' RAMSIN`
   sstdir=${sstdir%\'*}
   sstdir=${sstdir#*\'}
   
   # NDVIFILES
   ndvidir=`egrep '^   NDVIFPFX' RAMSIN`
   ndvidir=${ndvidir%\'*}
   ndvidir=${ndvidir#*\'}
   
   # ANALFILES
   analdir=`egrep '^   AFILEPREF' RAMSIN`
   analdir=${analdir%\'*}
   analdir=${analdir#*\'}
fi

##############################################################
# Execute repacking of HDF5 files

if [ -z "$3" ]; then
   if [ "$runtype" == "surface" ]; then
   
      # topo files
      lastmod=0 
      thisdir=$topdir
      repack
      
      # sfc files
      lastmod=0 
      thisdir=$sfcdir
      repack
   
      # sst files
      lastmod=0 
      thisdir=$sstdir
      repack
   
      # ndvi files
      lastmod=0 
      thisdir=$ndvidir
      repack
   
   elif [ "$runtype" == "analysis" ]; then
      
      # Copy RAMSIN to the directory on frost
      #scp RAMSIN ldgrant@frost.atmos.colostate.edu:$supcelldir 
#      scp RAMSIN $copydir


      # analysis files
      thisdir=$analdir
   
      # Determine job string for repack execution
      if [ $exectype -eq 1 ]; then         # LSF
         jobout=$(bjobs | grep $ramsname | grep RUN 2>&1 | sed "s/ .*//")
      elif [ $exectype -eq 2 ]; then       # Standard
         jobout=$(pidof -s $ramsname)
      fi
      echo "jobout1: $jobout"

      # Initialize filecounts.
      # The repacking process can take a long time. During the repacking, files
      # may be generated and the loop in 'repack' does not see them. These vars
      # will test if new files have been created.
      prefilecount=`ls -1 $thisdir*'.h5' | wc -l`
      postfilecount=0
      echo "prefilecount, postfilecount: $prefilecount $postfilecount"

      # Continue to execute repack if (1) the model is running OR 
      #                               (2) the number of files have changed since
      #                                   last repack, i.e. new files were generated
      lastmod=0
      while [[ -n "$jobout" ]] || [[ $prefilecount -ne $postfilecount ]]
      do   

         # Reset prefilecount before repack call
         prefilecount=`ls -1 $thisdir*'.h5' | wc -l`

         # Execute the repack function.
         repack

         # Regrab file count. If this number is different than prefilecount, the 
         # repack function will execute again
         postfilecount=`ls -1 $thisdir*'.h5' | wc -l`         
         echo "prefilecount, postfilecount: $prefilecount $postfilecount"

         # execute time controller
         sleep $checkint

         # Check again for job name to see if it finished
         # Determine job string for repack execution
         if [ $exectype -eq 1 ]; then         # LSF
            jobout=$(bjobs | grep $ramsname | grep RUN 2>&1)
         elif [ $exectype -eq 2 ]; then       # Standard
            jobout=$(pidof -s $ramsname)
         fi
         echo "jobout2: $jobout"
      done
      # One final check after while loop to catch the ending file
      #repack

      # We are done
      echo $ramsname" finished and x.repack.sh is terminating"
   else
      echo "Runtype not recognized!"
      echo "Runtype must be either 'surface' or 'analysis'"
      exit
   fi
else
   lastmod=0
   repacksecs=31556900 # seconds in a year. If files are older, change.
   thisdir=$3          # set directory to third argumemnt
   repack
   exit
fi

exit

###############################################################


