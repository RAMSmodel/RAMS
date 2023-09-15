# NOTES:
# For doing parallel simulations, compile MPICH or OPENMPI with both FORTRAN
# and C bindings turned on like the following configure statement for
# compiling MPICH. Most of the time you enter the software package main
# directory, issue the "configure" command, run "make", and then "make install".
# For the configure statement to work, your ".bash_profile" needs to know "which"
# C and Fortran compilers you are using. After compiling MPICH in this example,
# it will create "mpif90" and "mpicc" that you can use in this include.mk file
# for "F_COMP" and "C_COMP" env variables below. Make sure the paths to these
# compilers and associated libraries are in your .bash_profile PATH and
# LD_LIBRARY_PATH as needed.
# ./configure --prefix=/file-path-here-for-mpi/mpich-3.3.2 \
#  --enable-fast=O3 CPPFLAGS=-DNDEBUG CC=gcc FC=ifort

#############################################################################
# Define make (gnu make works best).
#############################################################################
MAKE=/usr/bin/make

#############################################################################
# Set your RAMS root path and version number.
#############################################################################
RAMS_ROOT=/home/smsaleeb/rams_git_dev
RAMS_VERSION=6.3.04

#############################################################################
# Set root locations for HDF5 I/O software.
# Choose parallel or serial processing option for your compile type.
# Typically can use "parallel" for either, but some supercomputers require
# use of the serial executable.
#############################################################################
HDF5_ROOT=/home/smsaleeb/software/hdf5-1.10.7

#############################################################################
# Set root locations for parallel processing MPI software.
# You can comment out MPI_ROOT for serial processing compile.
#############################################################################
MPI_ROOT=/home/smsaleeb/software/mpich-3.3.2

#############################################################################
# Do not change these 2. They point from RAMS_ROOT to the source code.
#############################################################################
MODEL=$(RAMS_ROOT)/src
UTILS_INCS=-I$(MODEL)/include

#############################################################################
# Set HDF libraries and include locations relative to HDF5_ROOT.
# Note that linking libraries below are in a particular order to work!
# You may not need to modify these 3 variables. Try this default.
# Everything after the "-Wl" identifies library location for shared objects.
# HDF5 requires libz (zip) and libsz (szip). These can be linked in LIBS
# in the compiler instead, but need to be on your computer system.
#############################################################################
#HDF5_LIBS=-L$(HDF5_ROOT)/lib -lhdf5_hl -lhdf5 \
#  -Wl,-rpath,/home/smsaleeb/software/szip-2.1/lib \
#  -Wl,-rpath,/home/smsaleeb/software/zlib-1.2.5/lib
HDF5_LIBS=-L$(HDF5_ROOT)/lib -lhdf5_hl -lhdf5
HDF5_INCS=-I$(HDF5_ROOT)/include
HDF5_DEFS=

#############################################################################
# TYPE OF COMPUTER SYSTEM (used for DEFINE statements for intrinsic calls).
# Do not change this variable.
#############################################################################
CMACH=PC_LINUX1  #Standard Linux (only option available now)

#############################################################################
# UNCOMMENT 1 LINUX FORTRAN COMPILER SET OF FLAGS BELOW
#############################################################################
# Note that some Fortran optimizations greater than -O1 will not necessarily
# produce identical results when running same executable multiple times.
# This is perhaps due to order of operations or unrolling loop order that
# changes solutions slightly. This is not ideal for research applications
# involving sensitivity studies that require exact duplication. As such,
# the Makefile is setup to allow different optimizations for different source
# files. The F_OPTS1 and F_OPTS2 variables below can hold different compiler
# flags as needed. Further, for duplication of results you may need to force
# IEEE standard which is done in examples below.
# Please coordinate variables below with the commands in the Makefile.
# Definitions of compiler flags below
# F_COMP      = path to compiler executable
# F_OPTS1     = Uses of flags for lower optimization for some files
# F_OPTS2     = Uses of flags for higher optimization for some files
# LOADER_OPTS = Flags and optimaztion for LOADER
# LIBS        = Extra system libraries for linking that some compilers need,
#   where the need for such libraries will be indicated by error messages
#   when compiling the code.

#*****************************
# FORTRAN INTEL IFORT COMPILER Single Precision
#*****************************
# (-g) for debugging, (-traceback) for more compiler error info
# (-check bounds) for array bounds checking, (-fp-model precise) for IEEE
# (-check uninit) for finding uninitialized variables, (-free) for free format
#F_COMP=/home/smsaleeb/intel/composer_xe_2011_sp1.8.273/bin/intel64/ifort
F_COMP=/home/smsaleeb/software/mpich-3.3.2/bin/mpif90
F_OPTS1=-free -O1 -fp-model precise
F_OPTS2=-free -O2 -fp-model precise
LOADER_OPTS= -free -O2 -fp-model precise
LIBS=-L/usr/lib/x86_64-linux-gnu -lrt -lpthread -lsz -lz

#*****************************
# FORTRAN INTEL IFORT COMPILER Double Precision
#*****************************
# (-g) for debugging, (-traceback) for more compiler error info
# (-check bounds) for array bounds checking, (-fp-model precise) for IEEE
# (-check uninit) for finding uninitialized variables, (-free) for free format
#F_COMP=/home/smsaleeb/intel/composer_xe_2011_sp1.8.273/bin/intel64/ifort
#F_OPTS1=-free -O1 -fp-model precise -real-size 64
#F_OPTS2=-free -O2 -fp-model precise -real-size 64
#LOADER_OPTS= -free -O2 -fp-model precise -real-size 64
#LIBS=-L/usr/lib/x86_64-linux-gnu -lrt -lpthread -lz -lsz

#*****************************
# FORTRAN PGI PGF90 COMPILER
#*****************************
# (-g) for debugging, (-Kieee) for IEEE, (-Mfree) for free format
#F_COMP=/opt/pgi/linux86-64/19.4/bin/pgf90
#F_OPTS1=-Mfree -O1 -Kieee
#F_OPTS2=-Mfree -O2 -Kieee
#LOADER_OPTS=-Mfree -O2 -Kieee
#LIBS=-L/usr/lib/x86_64-linux-gnu -lrt -lpthread -lz -lsz

#*****************************
# FORTRAN GFORTRAN COMPILER
# This compiler is not recommended and make not work as is. Gfortran seems
# to be less intuitive about certain code compared to PGF90 and IFORT.
#*****************************
# (-Wall) for warnings, (-ffree-form) for free format
# (-fno-sign-zero) for not making zeros negative values
# (-fcheck=bounds) check for array bounds issues
# (-fcheck=all) all runtime checking
#F_COMP=/usr/bin/gfortran
#F_OPTS1=-ffree-form -O1
#F_OPTS2=-ffree-form -O2
#LOADER_OPTS=-ffree-form -O2
#LIBS=-L/usr/lib/x86_64-linux-gnu -lrt -lpthread -lz -lsz

#############################################################################
# C compiler choice and flags (gcc) and (mpicc) are most common
# Add the "-DRAMS_DOUBLE_PREC" compiler flag to "C_OPTS" for turning on
# double precision rather than default single precision. Note that using
# double will make output and runtimes substantially longer but provide the
# extra precision needed in some highly sensitive simulations.
#
# Add the "-DENABLE_PARALLEL_COMPRESSION" flag to attempt doing parallel
# compression of HDF5 output. Sometime this works and sometimes it does not.
# But if it does work, it will make it past the first analysis file write.
# If it does not work, it will fail or hang on the first file write.
#
# Use -std=c99 if you need the c99 standard. Use -std=gnu99 potentially to
# stop unnecessary warnings related to "popen/pclose" in "dprep" code.
#
# "-w" turns off warnings, which in many cases are not a problem since our
# code has been well tested. New "gcc" versions sometime throw extra warnings
# that are not really as issue for us, but you can turn warnings back on by
# removing the "-w" if you wish to alter code to eliminate warnings.
#############################################################################
#C_COMP=gcc
C_COMP=/home/smsaleeb/software/mpich-3.3.2/bin/mpicc
C_OPTS=-O3 -DUNDERSCORE -DLITTLE -std=gnu99 -DENABLE_PARALLEL_COMPRESSION -w
#C_OPTS=-O3 -DUNDERSCORE -DLITTLE -std=gnu99 -DRAMS_DOUBLE_PREC \
#  -DENABLE_PARALLEL_COMPRESSION -w

#############################################################################
# System archive command syntax
# This might need different arguments depending on your Linux system. More
# recent versions have needed the "U" argument to prevent full RAMS recompile
# when only changing perhaps one file that does not have lots of dependencies
#############################################################################
ARCH=ar rsU

#############################################################################
# MPI - (Message Passing Interface) parallel processing choices required for
#  parallel enabled HDF5 and parallel processing on many-core systems.
#  Libraries to add to the PAR_LIBS list depend on your version of MPI
#  software and HDF5 software.
# Add or remove needed libraries to PAR_LIBS. Missing or needed libraries
#  will be highlighted in compile time error messages. Typical library
#  specs are: -lmpich, -lmpl, -lmpi, -lmpifort
# Comment out these "PAR_" lines for serial processing compile.
#############################################################################
PAR_INCS=-I$(MPI_ROOT)/include
#PAR_LIBS=-L$(MPI_ROOT)/lib -lmpich -lmpl
PAR_LIBS=-L$(MPI_ROOT)/lib
PAR_DEFS=-DRAMS_MPI
