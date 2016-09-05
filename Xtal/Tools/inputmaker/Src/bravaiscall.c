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
void bravaiscall(int index)	{
	
	float a,b,c;
	float alpha,beta,gamma;

	if ( index == 1|| index == 2 || index ==3 )	{
		a = get_a();
		b = a;
		c = a;
		alpha = 90.0000;
		beta = 90.0000;
		gamma = 90.0000;
	}
	else if ( index == 4 || index == 6 || index ==7 )	{
		a = get_a();
	        b = a;
		c = get_c();
                alpha = 90.0000;
		beta = 90.0000;
		gamma = 90.0000;

	}
	else if ( index == 5 ||index == -5 )	{
		a = get_a();
		b = a;
		c = get_c();
		alpha = get_alpha();
		beta = 90.0000;
		gamma = 90.0000;
	}
	else if ( index == 8 || index == 9 || index == -9 || index == 10 || index == 11)	{
		a = get_a();
		b = get_b();
		c = get_c();
		alpha = 90.0000;
		beta = 90.0000;
		gamma = 90.0000;
	}
	else if ( index == 12 || index == 13 )	{
		a = get_a();
		b = get_b();
		c = get_c();
		alpha = 90.0000;
		beta = 90.0000;	
		gamma = get_gamma();
	}
	else if ( index == -12 )	{
		a = get_a();
		b = get_b();
		c = get_c();
	 	alpha = 90.0000;
		beta = get_beta();
		gamma = 90.0000;
	}
	else if ( index == 14 )		{
		a = get_a();
		b = get_b();
		c = get_c();
		alpha = get_alpha();
		beta = get_beta();
		gamma = get_gamma();
	}
	else if ( index >14 )	{
		printf("\v");
		printf("\tError: Incorrect bravais index !\n");
		printf("\v");		
		exit(EXIT_SUCCESS);
	}
	writeoutput(index,a,b,c,alpha,beta,gamma);
}

