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

subroutine output( ipfile,opfile,bravais,a,b,c,alpha,beta,gamma,lat1,lat2,lat3, rlat1, rlat2, rlat3)
    !
    use kinds
    !
    use constants
    !
    !
    !
    character(len=32) :: ipfile,opfile
    integer:: bravais
    real(dp):: a,b,c
    real(dp):: alpha,beta,gamma
    real(dp),dimension(3):: lat1,lat2,lat3
    real(dp),dimension(3):: rlat1,rlat2,rlat3
    real(dp),dimension(3):: cross,rcross
    
    !
    open (unit=15,file=opfile)
        
    write(15,*)''
    write(15,'(a63)')'=--------------------------------------------------------='
    write(15,'(a63)')'=********************************************************='
    write(15,'(a63)')'=*                                                      *='
    write(15,'(a63)')'=*        XtalCalculator - Crys(X)tal Calculator        *='
    write(15,'(a63)')'=*               Copyright (C) 2013-2015                *='
    write(15,'(a63)')'=*                                                      *='
    write(15,'(a63)')'=********************************************************='
    write(15,'(a63)')'=--------------------------------------------------------='
    write(15,'(a63)')'=********************************************************='
    write(15,'(a63)')'=*                                                      *='
    write(15,'(a63)')'=*  Version      : 2.1.0-rev001                         *='
    write(15,'(a63)')'=*  Revision     : rev001                               *='
    write(15,'(a63)')'=*  Developer    : Anandavadivel Rajesh Prashanth       *='
    write(15,'(a63)')'=*  E-mail       : rajeshprasanth@rediffmail.com        *='
    write(15,'(a63)')'=*                                                      *='
    write(15,'(a63)')'=********************************************************='
    write(15,'(a63)')'=--------------------------------------------------------='
    write(15,*)''


    write(15,'(a19,a25)')'Input file  : ',ipfile
    write(15,'(a19,a25)')'Output file : ',opfile
        
    write(15,*)''
    write(15,'(a25)')'Data from input file'
    write(15,*)''

    if (bravais .eq. 1) then
        write(15,'(a40)')'Bravais lattice : Simple Cubic cell'

    else if (bravais .eq. 2) then
        write(15,'(a47)')'Bravais lattice : Face centered Cubic cell'

    else if (bravais .eq. 3) then
        write(15,'(a42)')'Bravais lattice : Body centered Cubic cell'

    else if (bravais .eq. 4) then
        write(15,'(a57)')'Bravais lattice : Simple Hexagonal and Trigonal cell'

    else if (bravais .eq. 5) then
        write(15,'(a52)')'Bravais lattice : Trigonal (3-fold axis c) cell'

    else if (bravais .eq. -5) then
        write(15,'(a56)')'Bravais lattice : Trigonal (3-fold axis <111>) cell'

    else if (bravais .eq. 6) then
        write(15,'(a45)')'Bravais lattice : Simple tetragonal cell'

    else if (bravais .eq. 7) then
        write(15,'(a52)')'Bravais lattice : Body centered tetragonal cell'

    else if (bravais .eq. 8) then
        write(15,'(a47)')'Bravais lattice : Simple Orthorhombic cell'

    else if (bravais .eq. 9) then
        write(15,'(a54)')'Bravais lattice : Base centered Orthorhombic cell'

    else if (bravais .eq. -9) then
        write(15,'(a54)')'Bravais lattice : Base centered Orthorhombic cell'

    else if (bravais .eq. 10) then
        write(15,'(a54)')'Bravais lattice : Face centered Orthorhombic cell'

    else if (bravais .eq. 11) then
        write(15,'(a54)')'Bravais lattice : Body centered Orthorhombic cell'

    else if (bravais .eq. 12) then
        write(15,'(a45)')'Bravais lattice : Simple Monoclinic cell'

    else if (bravais .eq. -12) then
        write(15,'(a45)')'Bravais lattice : Simple Monoclinic cell'

    else if (bravais .eq. 13) then
        write(15,'(a53)')'Bravais lattice : Body centered Monoclinic cell'

    else if (bravais .eq. 14) then
        write(15,'(a37)')'Bravais lattice : Triclinic cell'

    else
        write(*,'(a29)')'wrong bravais lattice!!!'
        stop
    end if

    call basic_vectors ( lat2, lat3, cross )
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
    vo = dot_product ( lat1, cross )    ! real cell volume
        !
    write(15,*)''
    write(15,'(a18)')'Intercepts : '
    write(15,*)''
    write(15,'(a14,f20.9,a6)')'a      : ',a,' (ang)'
    write(15,'(a14,f20.9,a6)')'b      : ',b,' (ang)'
    write(15,'(a14,f20.9,a6)')'c      : ',c,' (ang)'
    write(15,*)''
    write(15,'(a61)')'--------------------------------------------------------'
    write(15,*)''
    !
    !
    !
    write(15,'(a26)')'Interfacial angles : '
    write(15,*)''
    write(15,'(a13,f15.4,a6)')'alpha : ',alpha,' (deg)'
    write(15,'(a13,f15.4,a6)')'beta  : ',beta,' (deg)'
    write(15,'(a13,f15.4,a6)')'gamma : ',gamma,' (deg)'
    write(15,*)''
    write(15,'(a61)')'--------------------------------------------------------'
    write(15,*)''
        !
        !
        !
    write(15,'(a29)')'Primitive vectors (ang):'
    write(15,*)''
    write(15,'(3f20.9)')lat1
    write(15,'(3f20.9)')lat2
    write(15,'(3f20.9)')lat3
    write(15,*)''
    write(15,'(a23,f20.9,a8)')'Unit cell volume :',vo,'(ang^3)'
    write(15,*)''
    write(15,'(a61)')'--------------------------------------------------------'
    write(15,*)''
    !
    !
    !
    write(15,'(a30)')'Primitive vectors (bohr):'
    write(15,*)''
    write(15,'(3f20.9)')lat1*ang2bohr
    write(15,'(3f20.9)')lat2*ang2bohr
    write(15,'(3f20.9)')lat3*ang2bohr
    write(15,*)''
    write(15,'(a23,f20.9,a9)')'Unit cell volume :',vo*ang2bohr*ang2bohr*ang2bohr,' (bohr^3)'
    write(15,*)''
    write(15,'(a61)')'--------------------------------------------------------'
    write(15,*)''
    !
    !
    !
    call basic_vectors(rlat2,rlat3,rcross)

    rvol = dot_product(rlat1,rcross)

    write(15,'(a30)')'Reciprocal vectors (ang):'
    write(15,*)''
    write(15,'(3f20.9)')rlat1
    write(15,'(3f20.9)')rlat2
    write(15,'(3f20.9)')rlat3
    write(15,*)''
    write(15,'(a29,f20.9,a8)')'Reciprocal cell volume :',rvol,' (ang^3)'
    write(15,*)''
    write(15,'(a61)')'--------------------------------------------------------'
    !
    !
    !
    write(15,*)''
    write(15,'(a31)')'Reciprocal vectors (bohr):'
    write(15,*)''
    write(15,'(3f20.9)')rlat1*ang2bohr
    write(15,'(3f20.9)')rlat2*ang2bohr
    write(15,'(3f20.9)')rlat3*ang2bohr
    write(15,*)''
    write(15,'(a29,f20.9,a9)')'Reciprocal cell volume :',rvol*ang2bohr*ang2bohr*ang2bohr,' (bohr^3)'
    write(15,*)''
    write(15,'(a63)')'=--------------------------------------------------------='
    write(15,'(a29)')'Terminated sucessfully !'
    write(15,'(a63)')'=--------------------------------------------------------='


    return
end subroutine output
    

