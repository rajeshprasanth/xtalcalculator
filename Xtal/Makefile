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
#---------------------------------------------------------------------------
# This is the global makefile for XtalCalculator program.
# Do not modify this file.
#
#---------------------------------------------------------------------------
.SUFFIXES: .inc .f .f90
#
include make.inc

default :
	@(cd Scripts; bash splash.sh)
	@echo 
	@echo "         XtalCalculator Make targets "
	@echo 
	@echo "         Available targets"
	@echo
	@echo "         make xtalcalculator"
	@echo "         make resolve"
	@echo "         make tools"
	@echo "         make all"
	@echo "         make tar"
	@echo "         make clean"
	@echo "         make distclean"
	@echo "         make install"
	@echo "         make test"
	
	
splash:
	@(cd Scripts; bash splash.sh)
	
xtal :	
	@( cd Modules; $(MAKE) libmod.a )
	@( cd Src ; $(MAKE) xtalcalculator )
	@( cd Scripts ; bash compile.sh )
	@( cd Modules ; $(MAKE)  clean )
	@( cd Src ; $(MAKE) clean)
reso :
	@( cd Modules; $(MAKE) libmod.a )
	@( cd Src ; $(MAKE) resolve )
	@( cd Scripts ; bash compile.sh )
	@( cd Modules ; $(MAKE)  clean )
	@( cd Src ; $(MAKE) clean)	
	
tools:	
	@(cd Tools/inputmaker/Src; $(MAKE) inputmaker)
	#@( cd Tools/inputmaker/Src ; $(MAKE) clean)	

xtalcalculator : splash xtal

resolve : splash reso

all: splash xtal reso tools

tar:
	@( cd .. ; tar -zcvf XtalCalculator-2.1.0_rev001.tar.gz xtal )
	
clean :
	@( cd Modules ; $(MAKE) clean )		
	@( cd Src ; $(MAKE) clean)	
	@( cd Tools/inputmaker/Src ; $(MAKE) clean)	
	@( cd Bin ; /bin/rm -f *.o *.x *~ *.mod )
	/bin/rm -f *.o *.x *~ *.mod

distclean :
	@( cd Modules ; $(MAKE) clean )		
	@( cd Src ; $(MAKE) clean)	
	@( cd Tools/inputmaker/Src ; $(MAKE) clean)	
	/bin/rm -f *.o *.x *~ *.mod *.inc
	@( cd Bin ; /bin/rm -f *.o *.x *~ *.mod *.inc )	
	rmdir ./Bin

examples:



#install:


test:
	@(cd Testsuites; $(MAKE) all)
