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
.SUFFIXES: .h .c
#
include ../../../make.inc

INCLUDE= -I./

SRCS =  get_data.c	\
	bravaiscall.c	\
	writeoutput.c	\
	inputmaker.c	\
	
OBJS =  get_data.o	\
	bravaiscall.o	\
	writeoutput.o	\
	inputmaker.o	\
	
inputmaker : $(OBJS)
	$(CC) $(CFLAGS) -I./inputmaker.h $(OBJS) -o inputmaker.x

all: inputmaker

clean:
	-/bin/rm -f *.o *.x *~

#
#
include make.depend

