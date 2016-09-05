# -*- Autoconf -*-			                       	
#
# Copyright (C) 2013-2015 XtalCalculator (rajeshprasanth)
#
# This file is distributed under the terms of the GNU 
# General Public License. See the file `License' in the root 
# directory of the present distribution, or 
# http://www.gnu.org/copyleft/gpl.txt.
#
#------------------------------------------------------------------------------- 
#

#------------------------------#
#     Miscellaneous macros     #
#------------------------------#

AC_DEFUN([AC_SPLASH_DISPLAY],[
	echo " "
	echo "         =--------------------------------------------------------="
	echo "         =********************************************************="
	echo "         =*                                                      *="
	echo "         =*        XtalCalculator - Crys(X)tal Calculator        *="
	echo "         =*               Copyright (C) 2013-2015                *="
	echo "         =*                                                      *="
	echo "         =********************************************************="
	echo "         =--------------------------------------------------------="
	echo "         =********************************************************="
	echo "         =*                                                      *="
	echo "         =*  Version      : 2.1.0-rev001                         *="
	echo "         =*  Revision     : rev001                               *="
	echo "         =*  Developer    : Anandavadivel Rajesh Prashanth       *="
	echo "         =*  E-mail       : rajeshprasanth@rediffmail.com        *="
	echo "         =*                                                      *="
	echo "         =********************************************************="
	echo "         =--------------------------------------------------------="
	echo " "
])

AC_DEFUN([AC_MSG_DISPLAY],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  sec_title="$1"

  dnl Calculate amount of space chars needed for pretty-printing
  sec_spaces="                                                                      "
  sec_tmp="${sec_title}"
  while test "${sec_tmp}" != ""; do
    sec_spaces=`echo "${sec_spaces}" | sed -e 's/^.//'`
    sec_tmp=`echo "${sec_tmp}" | sed -e 's/^.//'`
  done

  echo ""
  echo " =============================================================================="
  echo " === ${sec_title}${sec_spaces} ==="
  echo " =============================================================================="
  echo ""
])
