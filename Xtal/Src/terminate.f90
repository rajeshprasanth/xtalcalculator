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

subroutine terminate()
    character(len=9) :: cdate, ctime

    call date_and_tim(cdate,ctime)

    write(*,'(a11,a9,a1,a9)')'End : ',cdate ,' ',ctime
    write(*,*)''
    write(*,'(a63)')'=--------------------------------------------------------='
    write(*,'(a29)')'Terminated sucessfully !'
    write(*,'(a63)')'=--------------------------------------------------------='
end subroutine terminate

