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

program POSCAR2XSF
	implicit none
	!
	!Variables
	!
	character(len=256) :: ipfile,opfile,option
	integer::error,i,status=0,j,k=1
        !
	!POSCAR Variables
	!
	character(len=128)::string,coordtype
	real::alat
	real::lattice(3,3)
	!
	logical::vasp5	
	integer,allocatable ::natoms(:)
	integer::atomtype=1,totalatoms=0
	character(len=2),allocatable:: atomsymbol(:)
	real,allocatable::atom_pos(:,:),d_atom_pos(:,:),c_atom_pos(:,:)
	!
 	if (iargc() .eq. 0) then
		write(*,*)''
	    	write(*,'(a41)')'POSCAR2XSF:: Fatal error : no input file'
        	write(*,'(a22)')'Conversion terminated'
		write(*,*)''
		stop
    	end if
	!	
	if (iargc() .eq. 2) then          !Using default VASP 5.X format
		call getarg(1,ipfile)
		call getarg(2,opfile)
		vasp5=.true.
	else	
		call getarg(1,option)
        	if (option .eq. "--vasp4") then
			vasp5=.false. 
			goto 111
		end if
		if (option .eq. "--help") then 
			call help() 
			stop
		else if (option .eq. "--version") then 
			call version() 
			stop
		else 
		     	write(*,*)'POSCAR2XSF ::Invalid option'
		     	call help() 
		     	stop
		end if
111	    	call getarg(2,ipfile)
	    	call getarg(3,opfile)	
    	end if
!--------------------------------------------------------------------------
	open(unit=5,file=trim(adjustl(ipfile)),form='formatted')
	rewind(5)
	read(5,'(a126)')string		!line 1 
	read(5,*)alat			!line 2
	read(5,*)lattice(1,:)		!line 3
	read(5,*)lattice(2,:)		!line 4
	read(5,*)lattice(3,:)		!line 5
!--------------------------------------------------------------------------		
	if (vasp5 .eqv. .true.) read(5,*)
	do while(status<=0)
		allocate(natoms(atomtype))
		read(5,*,iostat=status)(natoms(i),i=1,atomtype)
		atomtype=atomtype+1
		deallocate(natoms)
		backspace(5)
	end do
	atomtype = atomtype-2
  	allocate(atomsymbol(atomtype))
	allocate(natoms(atomtype))
	backspace(5)
	if (vasp5 .eqv. .true.) then
		backspace(5)
		read(5,*,iostat=status)(atomsymbol(i),i=1,atomtype)
	end if
      	
	read(5,*,iostat=status)(natoms(i),i=1,atomtype)
   
    	do i=1,atomtype
		totalatoms=natoms(i)+totalatoms
    	enddo

	read(5,*)coordtype		!line 8
	allocate(atom_pos(totalatoms,3))
!--------------------------------------------------------------------------
					
	if (coordtype(1:1)=='D' .or. coordtype(1:1)=='d') then
		do i=1,totalatoms
			read(5,*)atom_pos(i,:)
			d_atom_pos=atom_pos
			c_atom_pos=matmul(atom_pos,lattice)*alat
		end do
	end if
	if (coordtype(1:1)=='C' .or. coordtype(1:1)=='c' &
			.or. coordtype(1:1)=='K' .or. coordtype(1:1)=='k') then
		do i=1,totalatoms
			read(5,*)atom_pos(i,:)
			c_atom_pos=atom_pos
		end do
	end if
!---------------------------------------------------------------------------
!
!	do i=1,atomtype
!		do j=1,natoms(i)
!			write(*,*)adjustr(atomsymbol(i))
!		enddo
!	enddo
	open(unit=10,file=trim(adjustl(Opfile)),form='formatted')
	
	write(10,*)'# System name'
	write(10,*)'# ',adjustr(trim(string))
	write(10,*)'CRYSTAL'
	write(10,*)'PRIMVEC'
	write(10,'(3f15.8)')alat*0.529177*lattice(1,:)
	write(10,'(3f15.8)')alat*0.529177*lattice(2,:)
	write(10,'(3f15.8)')alat*0.529177*lattice(3,:)
	write(10,*)adjustr(trim('PRIMCOORD'))
	write(10,*)totalatoms,'1'
	do i=1,atomtype
		do j=1,natoms(i)
			write(10,'(a2,3f15.8)')adjustr(atomsymbol(i)),c_atom_pos(k,:)*0.529177
			k=k+1
		enddo
	enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!DO NOT TOUCH THIS PART! STILL UNDER DEVELOPMENT	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


	if (option .eq. "--verbose") then
		write(*,*)''
		end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program
!---------------------------------------------------------------------------------
subroutine help
	write(*,*)''
	write(*,'(a56)')'Usage: POSCAR2XSF [option] POSCAR_filename XSF_filename'
	write(*,'(a9)')'Option: '
	write(*,'(a49)')'--vasp4             Use VASP 4.X POSCAR format' 
	write(*,'(a42)')'(Default: VASP 5.X)'
	write(*,'(a47)')'--help              Display this information'
	write(*,'(a50)')'--version           Display version information'
	write(*,*)''
	write(*,*)'For reporting bugs, please contact:'
	write(*,*)'<rajeshprasanth@rediffmail.com>'
	write(*,*)''
	return    	
end subroutine

subroutine version
	write(*,*)''
	write(*,'(a23)')'POSCAR2XSF version 1.0'
	write(*,'(a50)')'Copyright (C) 2016 Rajesh Prashanth Anandavadivel'
	write(*,*)''
	write(*,'(a62)')'POSCAR2XSF is distributed in the hope that it will be useful,'
    	write(*,'(a63)')'but WITHOUT ANY WARRANTY; without even the implied warranty of'
    	write(*,'(a62)')'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
    	write(*,'(a45)')'GNU General Public License for more details.'
	write(*,*)''
	return
end subroutine

