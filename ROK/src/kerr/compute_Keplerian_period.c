/*==============================================================================
                                                                                       
   Given the orbital parameters associated with a GR Kerr geodesic (p,e)
                                                                                       
   this routine provides the associated Keplerian period, given by                     
                                                                                       
                              /    p    \  3/2                                         
                    T = 2 PI | --------- |                                             
                              \ 1 - e^2 /                                              
                                                                                       
-------------------------------------------------------------------------------
                                                                                       
   Last Update : 02.12.12 by PCM                                                              
==============================================================================*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "global_quantities_kerr.h"

int compute_Keplerian_period(double p_orbit, double eccentricity, double *t_kepler, long unsigned *n_one_t_kepler, double dt)
{
  double eccentricity2;                                                   // e^2
  double one_minus_ecc2;                                                  // 1 - e^2
  double p_over_one_minus_ecc2;                                           // p / ( 1 - e^2 ) 
  double p_over_one_minus_ecc2_to_the_3halfs;                             // [ p / ( 1 - e^2 ) ]^(3/2)


/*----- Some previous calculations -----*/
  eccentricity2 = eccentricity*eccentricity;                              // e^2
  one_minus_ecc2 = 1.0 - eccentricity2;                                   // 1 - e^2
  p_over_one_minus_ecc2 = p_orbit/one_minus_ecc2;                         // p / ( 1 - e^2 ) 
  p_over_one_minus_ecc2_to_the_3halfs = pow(p_over_one_minus_ecc2,1.5);   // [ p / ( 1 - e^2 ) ]^(3/2)


/*----- Keplerian period -----*/
  *t_kepler = 2.0*PI*p_over_one_minus_ecc2_to_the_3halfs;        
  
  
/*----- Number (ODE) Time steps within one Keplerian period -----*/  
  *n_one_t_kepler = (long unsigned) ceil( (*t_kepler)/dt );
  
  
/*----- Routine ends Here -----*/
  return 0;

}


