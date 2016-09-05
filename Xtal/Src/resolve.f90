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


!***********************************************************************************
!
!Main program for resolving the Primitive lattice vectors into the corresponding
!lattice parameters.
!
!***********************************************************************************

program resolve
    use kinds
    use constants

    character(len=4)::unit
    call display()

    write(*,'(a26)')'Type ang for Angstorm'
    write(*,'(a23)')'Type bohr for Bohr'
    write(*,'(a33)')'Enter the reference unit >>>'
    read(*,*)unit

    if (unit .eq. 'ang') then
        call convertang()

    else if (unit .eq. 'bohr') then
        call convertbohr()

    else
        write(*,'(a28)')'Error : Incorrect units'
        stop
    end if
    call terminate()
end program resolve

!***********************************************************************************
!
!Subroutine to resolve the lattice parameters if the inputs are in bohr.
!
!***********************************************************************************

subroutine convertbohr
            !
    use kinds
    use constants

    !
    !
    real(dp) :: alat
    real(dp) ::al,alph
    real(dp) ::bet,beta
    real(dp) ::gam,gamm
    real(dp),dimension(3) ::v1,v2,v3
    real(dp),dimension(3) ::ov1,ov2,ov3

    write(*,'(a24)')'Lattice constant = '
    read(*,*)alat

    write(*,'(a44)')'vector 1 (1),vector 1 (2),vector 1 (3):'
    read(*,*)v1(1),v1(2),v1(3)

    write(*,'(a44)')'vector 2 (1),vector 2 (2),vector 2 (3):'
    read(*,*)v2(1),v2(2),v2(3)

    write(*,'(a44)')'vector 3 (1),vector 3 (2),vector 3 (3):'
    read(*,*)v3(1),v3(2),v3(3)

    ov1(1) = v1(1)*alat
    ov1(2) = v1(2)*alat
    ov1(3) = v1(3)*alat

    ov2(1) = v2(1)*alat
    ov2(2) = v2(2)*alat
    ov2(3) = v2(3)*alat

    ov3(1) = v3(1)*alat
    ov3(2) = v3(2)*alat
    ov3(3) = v3(3)*alat

    write(*,'(a61)')'--------------------------------------------------------'
    write(*,*)''
    write(*,'(a22)')'Vectors(bohr) : '
    write(*,*)''
    write(*,'(3f20.9)')ov1
    write(*,'(3f20.9)')ov2
    write(*,'(3f20.9)')ov3
    write(*,*)''
    write(*,'(a61)')'--------------------------------------------------------'
    write(*,*)''

    write(*,'(a21)')'Vectors(ang) : '
    write(*,*)''
    write(*,'(3f20.9)')ov1*bohr2ang
    write(*,'(3f20.9)')ov2*bohr2ang
    write(*,'(3f20.9)')ov3*bohr2ang
    write(*,*)''
    write(*,'(a61)')'--------------------------------------------------------'
    write(*,*)''
    a=sqrt( (ov1(1)**2) + (ov1(2)**2) + (ov1(3)**2))
    b=sqrt( (ov2(1)**2) + (ov2(2)**2) + (ov2(3)**2))
    c=sqrt( (ov3(1)**2) + (ov3(2)**2) + (ov3(3)**2))

    write(*,'(a9,f20.9,a5)')'a = ',a,' bohr'
    write(*,'(a9,f20.9,a5)')'b = ',b,' bohr'
    write(*,'(a9,f20.9,a5)')'c = ',c,' bohr'
    write(*,*)''
    write(*,'(a61)')'--------------------------------------------------------'
    write(*,*)''
    write(*,'(a9,f20.9,a4)')'a = ',a*bohr2ang,' ang'
    write(*,'(a9,f20.9,a4)')'b = ',b*bohr2ang,' ang'
    write(*,'(a9,f20.9,a4)')'c = ',c*bohr2ang,' ang'
    write(*,*)''
    write(*,'(a61)')'--------------------------------------------------------'
    write(*,*)''

    al=( (dot_product(ov2,ov3)) / (b*c) )
    alph = acos(al)

    bet=( (dot_product(ov1,ov3)) / (a*c) )
    beta = acos(bet)

    gam=( (dot_product(ov1,ov2)) / (a*b) )
    gamm = acos(gam)

    write(*,'(a13,f20.9)')'alpha = ',(alph*180)/pi
    write(*,'(a13,f20.9)')'beta  = ',(beta*180)/pi
    write(*,'(a13,f20.9)')'gamma = ',(gamm*180)/pi
    write(*,*)''
    write(*,'(a61)')'--------------------------------------------------------'
    write(*,*)''

    write(*,'(a14,f20.9,a8)')'Volume = ',a*b*c*bohr2ang*bohr2ang*bohr2ang,' (ang^3)'
    write(*,'(a14,f20.9,a9)')'Volume = ',a*b*c,' (bohr^3)'
    write(*,*)''

    return
end subroutine convertbohr

!***********************************************************************************
!
!Subroutine to resolve the lattice parameters if the inputs are in angstrom.
!
!***********************************************************************************

subroutine convertang
            !
    use kinds
    use constants

    !
    !
    real(dp) :: alat
    real(dp) ::al,alph
    real(dp) ::bet,beta
    real(dp) ::gam,gamm
    real(dp),dimension(3) ::v1,v2,v3
    real(dp),dimension(3) ::ov1,ov2,ov3

    write(*,'(a24)')'Lattice constant = '
    read(*,*)alat

    write(*,'(a44)')'vector 1 (1),vector 1 (2),vector 1 (3):'
    read(*,*)v1(1),v1(2),v1(3)

    write(*,'(a44)')'vector 2 (1),vector 2 (2),vector 2 (3):'
    read(*,*)v2(1),v2(2),v2(3)

    write(*,'(a44)')'vector 3 (1),vector 3 (2),vector 3 (3):'
    read(*,*)v3(1),v3(2),v3(3)

    ov1(1) = v1(1)*alat
    ov1(2) = v1(2)*alat
    ov1(3) = v1(3)*alat

    ov2(1) = v2(1)*alat
    ov2(2) = v2(2)*alat
    ov2(3) = v2(3)*alat

    ov3(1) = v3(1)*alat
    ov3(2) = v3(2)*alat
    ov3(3) = v3(3)*alat

    write(*,'(a61)')'--------------------------------------------------------'
    write(*,*)''
    write(*,'(a21)')'Vectors(ang) : '
    write(*,*)''
    write(*,'(3f20.9)')ov1
    write(*,'(3f20.9)')ov2
    write(*,'(3f20.9)')ov3
    write(*,*)''
    write(*,'(a61)')'--------------------------------------------------------'
    write(*,*)''
    write(*,'(a22)')'Vectors(bohr) : '
    write(*,*)''
    write(*,'(3f20.9)')ov1*ang2bohr
    write(*,'(3f20.9)')ov2*ang2bohr
    write(*,'(3f20.9)')ov3*ang2bohr
    write(*,*)''
    write(*,'(a61)')'--------------------------------------------------------'
    write(*,*)''
    a=sqrt( (ov1(1)**2) + (ov1(2)**2) + (ov1(3)**2))
    b=sqrt( (ov2(1)**2) + (ov2(2)**2) + (ov2(3)**2))
    c=sqrt( (ov3(1)**2) + (ov3(2)**2) + (ov3(3)**2))

    write(*,'(a9,f20.9,a4)')'a = ',a,' ang'
    write(*,'(a9,f20.9,a4)')'b = ',b,' ang'
    write(*,'(a9,f20.9,a4)')'c = ',c,' ang'
    write(*,*)''
    write(*,'(a61)')'--------------------------------------------------------'
    write(*,*)''
    write(*,'(a9,f20.9,a5)')'a = ',a*ang2bohr,' bohr'
    write(*,'(a9,f20.9,a5)')'b = ',b*ang2bohr,' bohr'
    write(*,'(a9,f20.9,a5)')'c = ',c*ang2bohr,' bohr'
    write(*,*)''
    write(*,'(a61)')'--------------------------------------------------------'
    write(*,*)''

    al=( (dot_product(ov2,ov3)) / (b*c) )
    alph = acos(al)

    bet=( (dot_product(ov1,ov3)) / (a*c) )
    beta = acos(bet)

    gam=( (dot_product(ov1,ov2)) / (a*b) )
    gamm = acos(gam)

    write(*,'(a13,f20.9)')'alpha = ',(alph*180)/pi
    write(*,'(a13,f20.9)')'beta  = ',(beta*180)/pi
    write(*,'(a13,f20.9)')'gamma = ',(gamm*180)/pi
    write(*,*)''
    write(*,'(a61)')'--------------------------------------------------------'
    write(*,*)''

    write(*,'(a14,f20.9,a8)')'volume = ',a*b*c,' (ang^3)'
    write(*,'(a14,f20.9,a9)')'volume = ',a*b*c*ang2bohr*ang2bohr*ang2bohr,' (bohr^3)'
    write(*,*)''

    return
end subroutine convertang
