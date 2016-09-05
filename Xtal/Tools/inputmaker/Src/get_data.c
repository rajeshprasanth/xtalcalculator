/*!
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
*/
#include <stdio.h>
#include "inputmaker.h"
void get_tutorial(void)	{
	printf("\t+---------------------------------------------------+\n");	
	printf("\t|                     Tutorials                     |\n");
	printf("\t+-------------------------+-------------------------+\n");
	printf("\t|  1. Simple cubic cell   |  9. Base centered       |\n");
	printf("\t|                         |     orthorhombic cell   |\n");
	printf("\t+-------------------------+-------------------------+\n");	
	printf("\t|  2. Face centered       | -9. Base centered       |\n");
	printf("\t|     Cubic cell          |     orthorhombic cell-2 |\n");
	printf("\t+-------------------------+-------------------------+\n");	
	printf("\t|  3. Body centered       | 10. Face centered       |\n");
	printf("\t|     Cubic cell          |     orthorhombic cell   |\n");	
	printf("\t+-------------------------+-------------------------+\n");	
	printf("\t|  4. Simple hexagonal    | 11. Body centered       |\n");
	printf("\t|     and trigonal cell   |     orthorhomic cell    |\n");	
	printf("\t+-------------------------+-------------------------+\n");
	printf("\t|  5. Trigonal 3fold      | 12. Simple monoclinic   |\n");
	printf("\t|     axis c              |     unique axis c       |\n");	
	printf("\t+-------------------------+-------------------------+\n");	
	printf("\t| -5. Trigonal 3fold      |-12. Simple monoclinic   |\n");
	printf("\t|     axis <111>          |     unique axis b       |\n");	
	printf("\t+-------------------------+-------------------------+\n");	
	printf("\t|  6. Simple tetragonal   | 13. Monoclinic base     |\n");
	printf("\t|     cell                |     centered cell       |\n");	
	printf("\t+-------------------------+-------------------------+\n");	
	printf("\t|  7. Base centered       | 14. Triclinic cell      |\n");
	printf("\t|     tetragonal cell     |                         |\n");	
	printf("\t+-------------------------+-------------------------+\n");	
	printf("\t|  8. Simple orthogonal   |                         |\n");
	printf("\t|     cell                |                         |\n");
	printf("\t+---------------------------------------------------+\n");
}

int get_index(void )	{
	int index;
	printf("\tEnter bravais lattice number >>>\t");
	scanf("%d",&index);
	return (index);
}

float get_a(void)	{
	float a;	
	printf("\tEnter a (angstrom) >>>\t");
	scanf("%f",&a);
	return (a);
}

float get_b(void)	{
	float b;		
	printf("\tEnter b (angstrom) >>>\t");
	scanf("%f",&b);
	return (b);
}

float get_c(void)	{
	float c;
	printf("\tEnter c (angstrom) >>>\t");
	scanf("%f",&c);
	return (c);
}


float get_alpha(void)	{
	float alpha;
	printf("\tEnter alpha (degree) >>>\t");
	scanf("%f",&alpha);
	return (alpha);
}

float get_beta(void)	{
	float beta;
	printf("\tEnter beta (degree) >>>\t");
	scanf("%f",&beta);
	return (beta);
}

float get_gamma(void)	{
	float gamma;
	printf("\tEnter gamma (degree) >>>\t");
	scanf("%f",&gamma);
	return (gamma);
}








		
