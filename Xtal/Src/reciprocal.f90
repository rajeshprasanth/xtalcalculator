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

subroutine reciprocal( lat1, lat2, lat3, rlat1, rlat2, rlat3)
    !
    use kinds
    !
    use constants
    !
    !
    !
    real(dp),dimension(3)     ::lat1,lat2,lat3
    real(dp)                  :: rvol,a
    real(dp),dimension(3),intent(out)::rlat1,rlat2,rlat3
    real(dp)                  :: pre
    real(dp), dimension(3)    :: rcross
    real(dp), dimension(3)    :: cross,cross1,cross2,cross3
    !
    call basic_vectors (lat2,lat3,cross)
    call basic_vectors(lat2,lat3,cross1)
    call basic_vectors(lat3,lat1,cross2)
    call basic_vectors(lat1,lat2,cross3)

    a = dot_product (lat1,cross)
    pre = tpi/a

    rlat1 = pre *  cross1
    rlat2 = pre *  cross2
    rlat3 = pre *  cross3

    call basic_vectors(rlat2,rlat3,rcross)

    rvol = dot_product(rlat1,rcross)

    write(*,'(a30)')'Reciprocal vectors (ang):'
    write(*,*)''
    write(*,'(3f20.9)')rlat1
    write(*,'(3f20.9)')rlat2
    write(*,'(3f20.9)')rlat3
    write(*,*)''
    write(*,'(a29,f20.9,a8)')'Reciprocal cell volume :',rvol,' (ang^3)'
    write(*,*)''
    write(*,'(a61)')'--------------------------------------------------------'
    !
    !
    !
    write(*,*)''
    write(*,'(a31)')'Reciprocal vectors (bohr):'
    write(*,*)''
    write(*,'(3f20.9)')rlat1*ang2bohr
    write(*,'(3f20.9)')rlat2*ang2bohr
    write(*,'(3f20.9)')rlat3*ang2bohr
    write(*,*)''
    write(*,'(a29,f20.9,a9)')'Reciprocal cell volume :',rvol*ang2bohr*ang2bohr*ang2bohr,' (bohr^3)'
    write(*,*)''

    return
end subroutine reciprocal
