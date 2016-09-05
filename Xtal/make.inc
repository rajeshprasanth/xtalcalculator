# make.inc.  Generated from make.inc.in by configure.
#
# Copyright (C) 2013-2015 <Anandavadivel Rajesh Prashanth>
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.
#---------------------------------------------------------------------------
.SUFFIXES : .o .c .f90
#
#
#---------------------------------#
# .f.o and .c.o: do not modify    #
#---------------------------------#
#
.f90.o:
	$(FC) $(FCFLAGS) $(MODULES) -c $*.f90 -o $*.o
.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o
#

#
#---------------------------------#
# Flags for Fortran modules       #
#---------------------------------# 
#
MOD_FLAGS = -I
#
#---------------------------------#
# Fortran and C compilers         #
#---------------------------------#
#
FC = gfortran
CC = gcc
#
#---------------------------------#
# Fortran and C compilers Flags   #
#---------------------------------#
#
FCFLAGS = -g -O2 -ffree-form
CFLAGS = -g
#
#---------------------------------#
# Miscellaneous flags             #
#---------------------------------#
#
AR = ar
ARFLAGS = ruv
RANLIB = ranlib
INSTALL = /usr/bin/install -c
MAKEDIR = /bin/mkdir -p
LN_S = ln -s
SED = /bin/sed
AWK = gawk
GREP = /bin/grep
EGREP = /bin/grep -E
TAR = tar
WGET = wget
