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
subroutine basic_vectors ( vect1, vect2, outvect )
    !
    use kinds
    !
    implicit none
    real(dp),dimension(3):: vect1,vect2
    real(dp),dimension(3),intent(out):: outvect
    !
    outvect(1) = ((vect1(2)*vect2(3))-(vect2(2)*vect1(3)))
    outvect(2) = ((vect2(1)*vect1(3))-(vect1(1)*vect2(3)))
    outvect(3) = ((vect1(1)*vect2(2))-(vect2(1)*vect1(2)))
    !
    return
    stop
!
end subroutine basic_vectors
