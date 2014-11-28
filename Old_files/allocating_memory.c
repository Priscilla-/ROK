/*==============================================================================
                                                                                       
        THIS ROUTINE ALLOCATES MEMORY FOR THE DIFFERENT VARIABLES                            
                                                                                       
 -------------------------------------------------------------------------------
                                            Created by Priscilla on 20/12/2012
 -------------------------------------------------------------------------------
 Last Update : 28.12.12
===============================================================================*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h> 

#include "global_quantities_kerr.h"
#include "global_quantities_main.h"

int checks_errors(int);

int allocating_memory(void)
{
  
/*---- Allocating memory for storing the radial and polar coordinates on resonance ----*/
   // if ( (GV_Rres = (double *) malloc((size_t) ((GV_N_max)*sizeof(double)))) == NULL ) checks_main_errors(6);
    //if ( (GV_Thres = (double *) malloc((size_t) ((GV_N_max)*sizeof(double)))) == NULL ) checks_main_errors(6);
  
/*----- This is the end -----*/
  return 0;

}