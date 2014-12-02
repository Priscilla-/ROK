/*==============================================================================
 
  This routine modifies the computational parameters to fit ODE time steps
  within the sampling choosen, and also to fit the sampling with the total time 
  of evolution.
-------------------------------------------------------------------------------
                                            Created by Priscilla on 19/12/2012
-------------------------------------------------------------------------------
 Last Update : 09.07.13
==============================================================================*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "global_quantities_main.h"
#include "physical_quantities.h"
#include "global_quantities_kerr.h"

int modify_computational_parameters(void)
{
   
/*----- Setting Time Steps  -----*/

/**----- Nfreq  -----**/
  Nfreq = ceil( DTsec/(10.0*DT_INI*BHMASS_to_s) ); // Number of ODE Time Steps corresponding to the Sampling Frequency
 
/**----- Adjusting DT_INI to Nfreq -----**/
  DT_INI = 0.1*DTsec/BHMASS_to_s/((double) Nfreq);
		
/**----- NSAMPLING -----**/
  NSAMPLING = ceil(TIME_yrs*SECSYR/DTsec ) ;  // Total Number of Sampling Times along the inspiral

/**----- Adjusting TIME_yrs to NSAMPLING-----**/
  TIME_yrs = (NSAMPLING*DTsec)/SECSYR;  

/**----- NMAX -----**/   
  NMAX = Nfreq*NSAMPLING; //Total Number of ODE Time steps
	
  //printf("%lu %4.6e %lu %4.6e  %lu\n",Nfreq, DT_INI,  NSAMPLING, TIME_yrs, NMAX);exit(0);
    
/*----- This is the end -----*/
  return 0;

}




