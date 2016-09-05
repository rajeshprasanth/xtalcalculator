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
!  Description: A Subroutine for calculating lattice vectors for Trigonal 
!               (3-fold axis !!!<111>) cell.
!  Author:      Anandavadivel Rajesh Prashanth, rajeshprasanth@rediffmail.com
!___________________________________________________________________________________

subroutine trigonal_cell2 ( a, c, alpha, a1, a2, a3  )
    !
    use kinds
    !
    use constants
    !
    !=======================================================
    real(dp):: a,c,vo,tx,ty,tz,p,alpha,alphad
    real(dp):: ap,u,v
    real(dp),dimension(3),intent(out):: a1,a2,a3
    real(dp),dimension(3):: outvect,cross
    !
    !
    !
    alphad = ((pi*alpha)/180)
    p = cos (alphad)
    !
    ap = a/1.7320508075688772935274463415059
    tx = sqrt((1-p)/2)
    ty = sqrt((1-p)/6)
    tz = sqrt((1+(2*p))/3)
    !
    u = (tz)-(2*1.4142135623730950488016887242097*ty)
    v = (tz)+(1.4142135623730950488016887242097*ty)
    !
    a1(1) = ap*(u)
    a1(2) = ap*(v)
    a1(3) = ap*(v)
    !
    a2(1) = ap*(v)
    a2(2) = ap*(u)
    a2(3) = ap*(v)
    !
    a3(1) = ap*(v)
    a3(2) = ap*(v)
    a3(3) = ap*(u)
    !
    call basic_vectors ( a2, a3, cross )
    !================================================================
    !    volume = a1 * | a2 x a3 |
    !
    !         a1 = a1(1)  a1(2)  a1(3)
    !         a2 = a2(1)  a2(2)  a2(3)
    !         a3 = a3(1)  a3(2)  a3(3)
    !
    !            a1 = cellvector1
    !            a2 = cellvector2
    !            a3 = cellvector3
    !----------------------------------------------------------------
    vo = dot_product ( a1, cross )    ! real cell volume
    !
    write(*,*)''
    write(*,'(a25)')'Data from input file'
    write(*,*)''
    write(*,'(a56)')'Bravais lattice : Trigonal (3-fold axis <111>) cell'
    write(*,*)''
    write(*,'(a18)')'Intercepts : '
    write(*,*)''
    write(*,'(a14,f20.9,a6)')'a      : ',a,' (ang)'
    write(*,'(a14,f20.9,a6)')'b      : ',a,' (ang)'
    write(*,'(a14,f20.9,a6)')'c      : ',c,' (ang)'
    write(*,*)''
    write(*,'(a61)')'--------------------------------------------------------'
    write(*,*)''
    !
    !
    !
    write(*,'(a26)')'Interfacial angles : '
    write(*,*)''
    write(*,'(a13,f15.4,a6)')'alpha : ',alpha,' (deg)'
    write(*,'(a13,f15.4,a6)')'beta  : ',90.0000,' (deg)'
    write(*,'(a13,f15.4,a6)')'gamma : ',90.0000,' (deg)'
    write(*,*)''
    write(*,'(a61)')'--------------------------------------------------------'
    write(*,*)''
    !
    !
    !
    write(*,'(a29)')'Primitive vectors (ang):'
    write(*,*)''
    write(*,'(3f20.9)')a1
    write(*,'(3f20.9)')a2
    write(*,'(3f20.9)')a3
    write(*,*)''
    write(*,'(a23,f20.9,a8)')'Unit cell volume :',vo,'(ang^3)'
    write(*,*)''
    write(*,'(a61)')'--------------------------------------------------------'
    write(*,*)''
    !
    !
    !
    write(*,'(a30)')'Primitive vectors (bohr):'
    write(*,*)''
    write(*,'(3f20.9)')a1*ang2bohr
    write(*,'(3f20.9)')a2*ang2bohr
    write(*,'(3f20.9)')a3*ang2bohr
    write(*,*)''
    write(*,'(a23,f20.9,a9)')'Unit cell volume :',vo*ang2bohr*ang2bohr*ang2bohr,' (bohr^3)'
    write(*,*)''
    write(*,'(a61)')'--------------------------------------------------------'
    write(*,*)''
    return


end subroutine trigonal_cell2
