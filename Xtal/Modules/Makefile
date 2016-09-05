#
# Copyright (C) 2013-2015 <Anandavadivel Rajesh Prashanth>
#
# XtalCalculator is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.
#
# This file is distributed under the terms of the GNU General Public 
# License. See the file `License' in the root directory of the present 
# distribution, or http://www.gnu.org/copyleft/gpl.txt .
#---------------------------------------------------------------------------
.SUFFIXES: .inc .f .f90 
#
include ../make.inc

SRCS =  kinds.f90 \
	constants.f90 \
	global_variables.f90 \

OBJS = 	kinds.o \
 	constants.o \
	global_variables.o \

MODS = 	kinds.mod \
	constants.mod \
	global_variables.mod \

libmod.a :  $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	     $(RANLIB) $@

all :  $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	     $(RANLIB) $@

clean:
	-/bin/rm -f *.o *~ *.mod *.a
