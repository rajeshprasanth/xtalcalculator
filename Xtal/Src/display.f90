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

subroutine display
    !
    use constants
    !
    character(len=9) :: cdate, ctime

    call date_and_tim(cdate,ctime)
	
    write(*,*)''
    write(*,'(a63)')'=--------------------------------------------------------='
    write(*,'(a63)')'=********************************************************='
    write(*,'(a63)')'=*                                                      *='
    write(*,'(a63)')'=*        XtalCalculator - Crys(X)tal Calculator        *='
    write(*,'(a63)')'=*               Copyright (C) 2013-2015                *='
    write(*,'(a63)')'=*                                                      *='
    write(*,'(a63)')'=********************************************************='
    write(*,'(a63)')'=--------------------------------------------------------='
    write(*,'(a63)')'=********************************************************='
    write(*,'(a63)')'=*                                                      *='
    write(*,'(a63)')'=*  Version      : 2.2.0-rev001                         *='
    write(*,'(a63)')'=*  Revision     : rev001                               *='
    write(*,'(a63)')'=*  Developer    : Anandavadivel Rajesh Prashanth       *='
    write(*,'(a63)')'=*  E-mail       : rajeshprasanth@rediffmail.com        *='
    write(*,'(a63)')'=*                                                      *='
    write(*,'(a63)')'=********************************************************='
    write(*,'(a63)')'=--------------------------------------------------------='
    write(*,*)''

    write(*,'(a13,a9,a1,a9)')'Start : ',cdate ,' ',ctime
    write(*,*)''
    write(*,*)'                        (Cl-)-------------(Na+)'
    write(*,*)'                         /|                /|'
    write(*,*)'                        / |               / |'
    write(*,*)'                       /  |              /  |'
    write(*,*)'                    (Na+)-------------(Cl-) |'
    write(*,*)'                      |   |             |   |'
    write(*,*)'                      |   |     Xtal    |   |'
    write(*,*)'                      |   |  Calculator |   |'
    write(*,*)'                      |   |             |   |'
    write(*,*)'                      | (Na+)-----------|-(Cl-)'
    write(*,*)'                      |  /              |  /'
    write(*,*)'                      | /               | /'
    write(*,*)'                      |/                |/'
    write(*,*)'                    (Cl-)-------------(Na+)'
    write(*,*)''


    return
end subroutine display

