!
!    Copyright (C) <2013-2015>  <Anandavadivel Rajesh Prashanth>
!    
!    This file is part of XtalCalculator
!
!    XtalCalculator is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    XtalCalculator is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with XtalCalculator. If not, see <http://www.gnu.org/licenses/>.
!
!
!___________________________________________________________________________________
!
!  Author : Anandavadivel Rajesh Prashanth.
!  E-mail : rajeshprasanth@rediffmail.com
!___________________________________________________________________________________
 
subroutine argcheck(ipfile,opfile)
    !
    character(len=32) :: ifile,ofile
    character(len=32),intent(out) :: ipfile,opfile
    !
    integer::i

    if (iargc() .eq. 2) then
        call getarg(1,ifile)
        write(*,'(a19,a25)')'Input file  : ',ifile
        ipfile = ifile
        call getarg(2,ofile)
        write(*,'(a19,a25)')'Output file : ',ofile
        opfile = ofile
    else
        write(*,*)''    
        write(*,'(a63)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,'(a40)')'Error occured! Incorrect arguments!'
        write(*,'(a44)')'Usage :XtalCalculator file.in file.out!'
        write(*,'(a63)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,*)''
        stop
    end if

end subroutine argcheck
