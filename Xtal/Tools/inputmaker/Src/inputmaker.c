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
#include <string.h>
#include "inputmaker.h"
int main()	{

	long lind; 	
	char ch[8];
	int index;
	
	printf("\t=---------------------------------------------------=\n");	
	printf("\t                    INPUT MAKER                      \n");
	printf("\t=---------------------------------------------------=\n");
	printf("\v");
	do	{
		printf("\tDo you need tutorial ? (yes/no)>>>\t");
		scanf("%s",ch);
		printf("\v");
	
		if (strcmp(ch, "yes") == 0)	{
    			get_tutorial();
			lind = 1;
		}
		else	if(strcmp(ch, "no") == 0)	{
				index = get_index();
				bravaiscall(index);
		}		
		else	{
			printf("\v"); 
			printf("\n\tError: Incorrect option! Try again.\n");
			printf("\v");
			exit (EXIT_SUCCESS);
		}

	}while(lind>0);
 	return (EXIT_SUCCESS);
}

