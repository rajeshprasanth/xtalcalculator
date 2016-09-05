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

subroutine bravais_check(bravais,a,b,c,alpha,beta,gamma,a1,a2,a3)
    !
    use kinds
    !
    integer,intent(in):: bravais
    real(dp),intent(in):: a,b,c
    real(dp),intent(in):: alpha,beta,gamma
    real(dp),dimension(3),intent(out):: a1,a2,a3

    if (bravais .eq. 1 ) then
        call cubic_cell1 (a,a1,a2,a3)

    else    if (bravais .eq. 2) then
        call cubic_cell2 (a,a1,a2,a3)

    else    if (bravais .eq. 3) then
        call cubic_cell3 (a,a1,a2,a3)

    else    if (bravais .eq. 4) then
        call hextri_cell1 (a,c,a1,a2,a3)

    else    if (bravais .eq. 5) then
        call trigonal_cell1 ( a, c, alpha,a1, a2, a3  )

    else    if (bravais .eq. -5) then
        call trigonal_cell2 ( a, c, alpha,a1, a2, a3  )

    else    if (bravais .eq. 6) then
        call tetra_cell1 ( a, c, a1, a2, a3  )

    else    if (bravais .eq. 7) then
        call tetra_cell2 ( a, c, a1, a2, a3  )

    else    if (bravais .eq. 8) then
        call ortho_cell1 ( a, b, c, a1, a2, a3  )

    else    if (bravais .eq. 9) then
        call ortho_cell2 ( a, b, c, a1, a2, a3  )

    else    if (bravais .eq. -9) then
        call ortho_cell3 ( a, b, c, a1, a2, a3  )

    else    if (bravais .eq. 10) then
        call ortho_cell4 ( a, b, c, a1, a2, a3  )

    else    if (bravais .eq. 11) then
        call ortho_cell5 ( a, b, c, a1, a2, a3  )

    else    if (bravais .eq. 12) then
        call mono_cell1( a,b,c,gamma,a1,a2,a3 )

    else    if (bravais .eq. -12) then
        call mono_cell2( a,b,c,beta,a1,a2,a3 )

    else    if (bravais .eq. 13) then
        call mono_cell3( a,b,c,gamma,a1,a2,a3 )

    else    if (bravais .eq. 14) then
        call tri_cell( a,b,c,alpha,beta,gamma,a1,a2,a3 )

    else
        write(*,'(a29)')'wrong bravais lattice!!!'
        stop

    end if

end subroutine bravais_check
