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

module constants
    !
    use kinds, only : dp
    !
    implicit none
    !
    save
    !
    real(dp), parameter :: pi = 3.14159265358979323846_DP
    real(dp), parameter :: tpi = 2.0_DP * pi
    real(dp), parameter :: fpi = 4.0_DP * pi
    real(dp), parameter :: ang2bohr = 1.8897268777743552724324753343399_DP
    real(dp), parameter :: bohr2ang = 0.529177_DP
!
end module constants
