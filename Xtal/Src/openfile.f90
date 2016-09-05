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

subroutine openfile(ipfile,opfile)
    !
    integer ::ierror,oerror
    character(len=32) :: ipfile,opfile
    !
    open (unit=13,file=ipfile,status="old",iostat=ierror)

    if (ierror /= 0) then
        write(*,*)''
        write(*,'(a51)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,'(a41)')'Error occured! Input file not found!'
        write(*,'(a36)')'Check if Input file is present!'
        write(*,'(a52)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        stop
    end if

    open (unit=15,file=opfile,iostat=oerror)

    if (oerror /= 0) then
        write(*,*)''
        write(*,'(a51)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,'(a49)')'Error occured! Unable to create Output file!'
        write(*,'(a45)')'Check if sufficient memory is available!'
        write(*,'(a29)')'Delete if already exist!'
        write(*,'(a51)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        stop
    end if

end subroutine openfile
