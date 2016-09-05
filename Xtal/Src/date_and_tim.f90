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
    
subroutine date_and_tim (cdate, ctime)
    !
    implicit none
    character (len=9) :: cdate, ctime
    !
    character(len=3), dimension(12) :: months
    data months /'Jan','Feb','Mar','Apr','May','Jun',                     &
        'Jul','Aug','Sep','Oct','Nov','Dec'/
    INTEGER date_time(8)
    !
    call date_and_time(values=date_time)
    !
    write (cdate,'(i2,a3,i4)') date_time(3), months(date_time(2)), date_time(1)
    write (ctime,'(i2,":",i2,":",i2)') date_time(5), date_time(6), date_time(7)

end subroutine date_and_tim
