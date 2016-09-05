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

MODULES = $(MOD_FLAGS)../Modules

MODLIB = ../Modules/libmod.a


SRCS1 = date_and_tim.f90 \
	display.f90 \
	argcheck.f90 \
	openfile.f90 \
	readfile.f90 \
	bravais_check.f90 \
	basic_vector.f90 \
	cubic_cell1.f90	\
	cubic_cell2.f90	\
	cubic_cell3.f90	\
	hextri_cell1.f90 \
	trigonal_cell1.f90 \
	trigonal_cell2.f90 \
	tetra_cell1.f90 \
	tetra_cell2.f90 \
	ortho_cell1.f90 \
	ortho_cell2.f90 \
	ortho_cell3.f90 \
	ortho_cell4.f90 \
	ortho_cell5.f90 \
	mono_cell1.f90 \
	mono_cell2.f90 \
	mono_cell3.f90 \
	tri_cell.f90 \
	reciprocal.f90 \
	output.f90 \
	terminate.f90 \
	xtalcalculator.f90 \
	
	  
OBJS1 = date_and_tim.o	\
	display.o	\
	argcheck.o	\
	openfile.o	\
	readfile.o	\
	bravais_check.o	\
	basic_vector.o	\
	cubic_cell1.o	\
	cubic_cell2.o	\
	cubic_cell3.o	\
	hextri_cell1.o	\
	trigonal_cell1.o \
	trigonal_cell2.o \
	tetra_cell1.o \
	tetra_cell2.o \
	ortho_cell1.o \
	ortho_cell2.o	\
	ortho_cell3.o \
	ortho_cell4.o \
	ortho_cell5.o \
	mono_cell1.o \
	mono_cell2.o \
	mono_cell3.o \
	tri_cell.o \
	reciprocal.o	\
	output.f90 \
	terminate.o \
	xtalcalculator.o \
	
	
SRCS2=  display.f90	\
	date_and_tim.f90 \
	terminate.f90 \
	resolve.f90 \

OBJS2=  display.o \
	date_and_tim.o \
	terminate.o \
	resolve.o \

xtalcalculator : $(OBJS1)
	$(FC) $(FCFLAGS) $(OBJS1) $(MODULES) $(MODLIB)  -o xtalcalculator.x

resolve : $(OBJS2)
	$(FC) $(FCFLAGS) $(OBJS2) $(MODULES) $(MODLIB) -o resolve.x

all: $(OBJS1) $(OBJS2)  
	$(FC) $(FCFLAGS) $(OBJS1) $(MODULES) $(MODLIB)  -o xtalcalculator.x
	$(FC) $(FCFLAGS) $(OBJS2) $(MODULES) $(MODLIB) -o resolve.x
	 
clean:
	-/bin/rm -f *.o *.x *~ *.mod
