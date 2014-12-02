/*==============================================================================
 
    Clossing file: clean up and freen memory
 -------------------------------------------------------------------------------
                                            Created by Priscilla on 19/12/2012
 -------------------------------------------------------------------------------
 Last Update : 14.01.13
 =============================================================================*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<sys/types.h>
#include<time.h>

#include "global_quantities_kerr.h"
#include "global_quantities_main.h"

void free_Realmatrix(double **,const long , const long,const long , const long );
void free_Realvector(double *v, const long nl, const long nh);

int wrap_up()
{
    
 
/*----- Freeing Memory -----*/
  free_Realmatrix(E_p, 1, 4, 1, 4);
  free_Realmatrix(E_x, 1, 4, 1, 4);
  
   
  
  printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  printf("                      =============                   \n");
  printf("                      RUN  FINISHED                   \n");
  printf("                      =============                  \n\n");
  printf(" http://grooveshark.com/s/This+Is+The+End/4TXizK?src=5 \n\n");
  printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  printf("\n");
  printf("\n");
  printf("\n");
/*----- This is the end -----*/
  return 0;
}

