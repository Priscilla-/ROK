/*==============================================================================
 Given the current values of the constants of motion, (E,Lz,C/Q), the values of
 their evolution, (E_dot,Lz_dot,C_dot), and the value of the period of coordinate 
 time that has passed [i.e. the total
 coordinate time for which we have evolved the Geodesic equations with the 
 (E,Lz,C/Q) fixed], computes the new constants of motion.
  NOTE: See below for further information
-------------------------------------------------------------------------------                                                                                       
 Last Update : 30.06.13  by PCM                                                            
===============================================================================*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "global_quantities_kerr.h"


int compute_new_constants_of_motion(double time_elapsed,
                                    double energy_dot, double angular_momentum_z_dot, double q_carter_dot,
                                    double *energy, double *angular_momentum_z, double *c_carter, double *q_carter)
{
  double dt_energy, dt_angular_momentum_z, dt_q_carter;
 
/*----- Computing the changes in the constants of motion (E,Lz,C) -----*/
  dt_energy = energy_dot*time_elapsed;
  
  dt_angular_momentum_z = angular_momentum_z_dot*time_elapsed;
  
  dt_q_carter = q_carter_dot*time_elapsed;
  
 
/*----- Compute the new constants of motion (E_new,Lz_new,C_new/Q_new) -----*/
  *energy += dt_energy;
  
  *angular_momentum_z += dt_angular_momentum_z;
  
  *q_carter += dt_q_carter;
  
  *c_carter = (*q_carter) + ( (*angular_momentum_z) - SPIN*(*energy) )*( (*angular_momentum_z) - SPIN*(*energy) );
    

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
/*   Q = C - ( a E - Lz )^2                                                              */
/*                                                                                       */
/*****************************************************************************************/
/*                                                                                       */
/*   Computing the new constants of motion using the following simple procedure:         */
/*                                                                                       */
/*    1. Compute the changes in the constants of motion (E,Lz,C):                        */
/*                                                                                       */
/*        DE  = E_dot x T        ,                                                       */
/*                                                                                       */
/*        DLz = Lz_dot x T        ,                                                      */
/*                                                                                       */
/*        DQ  = Q_dot x T        .                                                       */
/*                                                                                       */
/*                                                                                       */
/*    2. Compute the new constants of motion (E_new,Lz_new,C_new/Q_new):                 */
/*                                                                                       */
/*        E_new  = E_current + DE ,                                                      */
/*                                                                                       */
/*        Lz_new = Lz_current + DLz ,                                                    */
/*                                                                                       */
/*        C_new  = C_current + DC ,                                                      */
/*                                                                                       */
/*        Q_new  = C_new - (Lz_new - a E_new)^2 .                                        */
/*                                                                                       */
/*   Details of the computation of (E_dot,Lz_dot,C_dot) can be found in the routine:     */
/*                                                                                       */
/*        compute_ELCQ_dot.c                                                             */
/*                                                                                       */
/*****************************************************************************************/
