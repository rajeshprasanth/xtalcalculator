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
!Main driver program for XtalCalculator.
!
!***********************************************************************************

program xtalcalculator

    use kinds
    use global_variables

    !***********************************************************************************
    !
    !Prints Program name, version, date and time of start of run and a beautiful crystal
    !lattice of Sodium Chloride.
    !
    !***********************************************************************************
    
    call display()


    !***********************************************************************************
    !
    !Check the command line arguments of XtalCalculator and returns input and
    !output file name. If the arguments are incorrect then returns error
    !message.
    !
    !***********************************************************************************

    call argcheck(ipfile,opfile)

    !***********************************************************************************
    !
    !Opens the input and output files. If either of files is found to be
    !missing then returns error message.
    !
    !***********************************************************************************

    call openfile(ipfile,opfile)

    !***********************************************************************************
    !
    !Read the input file content and returns the bravais lattice number,a(ang),b(ang),c(ang)
    !alpha(deg),beta(deg) and gamma(deg).
    !
    !***********************************************************************************

    call readfile(ipfile,opfile,bravais,a,b,c,alpha,beta,gamma)

    !***********************************************************************************
    !
    !Resolve bravais lattice number to bravais lattice and call corresponding bravais
    !lattice subroutine for computing primitive cell vectors (in ang and bohr) and cell
    !volume (in ang^3 and bohr^3).
    !
    !***********************************************************************************

    call bravais_check(bravais,a,b,c,alpha,beta,gamma,a1,a2,a3)

    !***********************************************************************************
    !
    !Calculate reciprocal lattice vectors (in ang and bohr) and reciprocal cell volume or
    !Brillouin zone volume (in ang^3 and bohr^3) from the primitive lattice vectors.
    !
    !***********************************************************************************

    call reciprocal( a1, a2, a3, b1, b2, b3)

    !***********************************************************************************
    !
    !Print all the computed values in the output file with the output file name specified
    !in the command-line arguments of Xtalcalculator.
    !
    !***********************************************************************************

    call output( ipfile,opfile,bravais,a,b,c,alpha,beta,gamma,a1, a2, a3, b1, b2, b3)

    !***********************************************************************************
    !
    !Prints the date and time of end of run and termination message.
    !
    !***********************************************************************************

    call terminate()

!***********************************************************************************
!
!End of program
!
!***********************************************************************************

end program xtalcalculator
