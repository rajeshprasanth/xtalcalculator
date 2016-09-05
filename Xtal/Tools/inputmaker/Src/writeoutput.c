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
#include <stdlib.h>
#include "inputmaker.h"
void writeoutput(int index,float a,float b,float c,float alpha,float beta,float gamma)	{
	
	FILE *fp;
	char filename[256];
	
	printf("\v");
	printf("\tEnter output file name >>>\t");
	scanf("%s",filename);
		
	printf("\v");
		
	fp = fopen(filename,"w");
	if (fp == NULL)	{
		printf("\tFile cannot be created!\n");
		exit (EXIT_SUCCESS);
	}
	
	fprintf(fp, "\t %d \t                 !\tindex bravais lattice index\n",index);
	fprintf(fp, "\t a\t %f \t !\ta (angstrom)\n",a);
	fprintf(fp, "\t b\t %f \t !\tb (angstrom)\n",b);
	fprintf(fp, "\t c\t %f \t !\tc (angstrom)\n",c);

	fprintf(fp, "\t alpha\t %f \t !\talpha  (degree)\n",alpha);
	fprintf(fp, "\t beta\t %f  \t !\tbeta   (degree)\n",beta);
	fprintf(fp, "\t gamma\t %f \t !\tgamma  (degree)\n",gamma);

	fclose(fp);
}
