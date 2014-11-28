/*==============================================================================
   Given the new values of the constants of motion, (E_new,Lz_new,C_new/Q_new),
   this function checks whether the new Geodesic Orbit is a plunge Orbit or
   it is a regular Geodesic.
 -------------------------------------------------------------------------------
  Last Update : 29.06.13    by PCM       
==============================================================================*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "global_quantities_kerr.h"
int checks_errors(int);

//int check_if_orbit_plunges(char *FLAG_iteration, double eccentricity, double q_carter)
int check_if_orbit_plunges(double eccentricity, double q_carter)

{
     
/*----- Checking whether the Carter constant has become negative -----*/
  if ( q_carter < 0.0 )
  {
      checks_errors(20);
  }//*FLAG_iteration = '0';
  
  
/*----- Checking whether the eccentricity has become negative -----*/
  if ( eccentricity < 0.0 )
  {
      checks_errors(20);
  }//  *FLAG_iteration = '0';
	
/*----- Checking whether we have reached separatrix-----*/
  if(R_peri == R_3)
  {
        printf("[!] We have hit the separatrix"); exit(0);
  }
  
/*----- Routine ends Here -----*/
  return 0;

}

/*****************************************************************************************/
/*                                                                                       */
/*   The conventions of the ROK Code for the Carter Constant are:                        */
/*                                                                                       */
/*      4 .2       2   2           2                         2   2                       */
/*   rho  r  = [ (r + a )E - a Lz ]  - Delta [ Q + (a E - Lz) + r  ]                     */
/*                                                                                       */
/*   Q = C + ( a E - Lz )^2                                                              */
/*                                                                                       */
/*****************************************************************************************/
