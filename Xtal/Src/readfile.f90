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

subroutine readfile(ipfile,opfile,bravais,a,b,c,alpha,beta,gamma)

    use kinds
    character(len=32) :: ipfile,opfile

    character(len=2)::label_a,label_b,label_c
    character(len=5)::label_alpha,label_beta,label_gamma


    integer,intent(out):: bravais
    real(dp),intent(out):: a,b,c
    real(dp),intent(out):: alpha,beta,gamma

    open (unit=13,file=ipfile)


    read(13,*)bravais
    !
    !
    !
    !
    !
    read(13,*)label_a,a
    if ( trim(label_a) /= 'a')  then
        write(*,*)''
        write(*,'(a51)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,'(a39)')'Error  : Incorrect label in line 2'
        write(*,'(a51)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,*)''
        stop
    end if
    !
    read(13,*)label_b,b
    if ( trim(label_b) /='b')  then
        write(*,*)''
        write(*,'(a51)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,'(a39)')'Error  : Incorrect label in line 3'
        write(*,'(a51)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,*)''
        stop
    end if
    !
    read(13,*)label_c,c
    if ( trim(label_c) /='c')  then
        write(*,*)''
        write(*,'(a51)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,'(a39)')'Error  : Incorrect label in line 4'
        write(*,'(a51)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,*)''
        stop
    end if
    !
    !
    !
    !
    !
    read(13,*)label_alpha,alpha
    write(*,*)
    if ( trim(label_alpha) /='alpha')  then
        write(*,*)''
        write(*,'(a51)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,'(a39)')'Error  : Incorrect label in line 5'
        write(*,'(a51)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,*)''
        stop
    end if
    !
    read(13,*)label_beta,beta
    if ( trim(label_beta) /='beta')  then
        write(*,*)''
        write(*,'(a51)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,'(a39)')'Error  : Incorrect label in line 6'
        write(*,'(a51)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,*)''
        stop
    end if
    !
    read(13,*)label_gamma,gamma
    if ( trim(label_gamma) /='gamma')  then
        write(*,*)''
        write(*,'(a51)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,'(a39)')'Error  : Incorrect label in line 7'
        write(*,'(a51)')'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,*)''
        stop
    end if
    close(13)
end subroutine readfile
