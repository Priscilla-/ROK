/*==============================================================================
 
 This routine evolves the geodesic equations for the Kerr Spacetime using the 
 system of equations that one derives when using the first integrals associated
 with the timelike and spatial axial Killing vectors and the one associated 
 with the Killing tensor (which makes the system to be integrable).  
 To that end, Boyer-Lindquist coordinates are used and, to avoid the numerical 
 problems associated with their turning points, we introduce two more angular 
 coordinates, chi and psi, associated with theta and r, respectively.  
                                                                                         
  This Routine:
  @ Computes the fundamental frequencies associated with the SCO geodesic

  @ Evolves the equations for (T_p,Psi_p,Chi_p,Phi_p) and saves their
    values in a file.
 
  NOTE: See below for further information
 -------------------------------------------------------------------------------
                                            Created by Priscilla on 20/12/2012
 -------------------------------------------------------------------------------
 Last Update: 22.09.13
 =============================================================================*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/types.h>
#include<time.h>

#include "global_quantities_kerr.h"
#include "global_quantities_main.h"
#include "physical_quantities.h"
#include "macros_2.h"

int check_if_orbit_plunges(double, double);
int checks_errors(int);
int compute_ELCQ_dot(double, double, double, double, double, double, double, double, double *, double *, double *);
int compute_Keplerian_period(double, double, double *, long unsigned *);
int compute_new_constants_of_motion(double, double, double, double, double *, double *, double *, double *);
int ELCQ_to_peiz(double, double, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int geodesic_evolution(double *, double *, double *,double *, double, double, double, double, double);
int pez_to_ELCQ(char, double, double, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *);

int evolution_RR(FILE *data )
{
  char FLAG_circular;                            // This flags tells whether the orbit has to be considered circular or not
  char FLAG_iteration;                           // Flag that indicates whether Radiation Reaction Effects have driven the orbit to a plunge trajectory
 
  long unsigned n_global;                        // Global Time Step Counter 
  long unsigned n_kepler_max;                    // Total Number of Time Steps within a Time Interval without applying Radiation Reaction Effects
  long unsigned n_sampling;                      // Counter of Sampling Time

  double t_initial;                              // Initial Time for ODE Integrator
  double t_in_BHmass, t_in_sec;                  // evolution times
  double time_elapsed;                           // idem
  double r_peri_0e,r_apo_0e;                     // Quantities for the particular case eccentricity = 0 [Useful to correct the Fluxes near eccentricity = 0]
  double inclination_0e;                         // idem
  double theta_minus_0e,z_minus_0e;              // idem
  double energy_0e,angular_momentum_z_0e;        // idem
  double c_carter_0e,q_carter_0e;                // idem
  double lambdat;                                //Mino time parameterized in funcion of the coordinate time
    
/*----- Initializing  Variables and Memory Allocation -----*/
  FLAG_iteration = '1';                 // Initializing the Flag about the Radiation Reaction Iteration Process
  n_sampling = 0;                       // Initialiazing the Sampling Time Counter
  n_kepler_max = N_periods*StepKT;      // Total Number of Time Steps within a Time Interval without applying Radiation Reaction Effects
    
  lambdat = 0.0;
  
   
/*======= EVOLUTION  LOOP  STARTS  HERE! ========*/
  for (n_global = 0; n_global <= NMAX; n_global++)
  {

/*----- Adjusting  times -----*/
    t_in_BHmass = n_global*DT_MBH; 
    t_in_sec    = t_in_BHmass*BHMASS_to_s;
      
/*----- Computing the Time Elapsed between applying  RR -----*/
    time_elapsed += DT_MBH;
      
/*----- Setting the Time, Time Intervals, and Variables for ODE Evolution of the SCO trajectory -----*/		
    t_initial = t_in_sec/BHMASS_to_s - DT_MBH;
	  
    CHI_NEW = CHI_OLD;
    PSI_NEW = PSI_OLD;
    PHI_NEW = PHI_OLD;
     
/*---- Computing radial and polar Boyer-Lindquist's coordinates ----*/
    R_p_o = P_p/(1.0 +ECCENTRICITY*cos(PSI_NEW) );
    COSTH_p_o = cos(THETA_min)*cos(CHI_NEW);

       
/*----- ODE Geodesic Evolution -----*/
    geodesic_evolution(&PSI_NEW, &CHI_NEW, &PHI_NEW,&lambdat, t_initial, t_in_BHmass, DT_INI, DT_MIN, ACCURACY);
    
/*----- Radiation Reaction Computations [UPDATING THE ORBITAL PARAMETERS/CONSTANTS OF MOTION] -----*/
  if ( (FLAG_RR == 'y')&& (n_global%n_kepler_max == 0) )
  {
/**----- Compute the Constants of Motion and related Quantities for Eccentricity = 0 -----**/
//	  FLAG_circular = 'y';
    pez_to_ELCQ(FLAG_circular, P_p, 0.0, THETA_inc, SPIN,
                &inclination_0e, &r_peri_0e, &r_apo_0e, &theta_minus_0e, &z_minus_0e,
                &energy_0e, &angular_momentum_z_0e, &c_carter_0e, &q_carter_0e);

/**----- Compute the Time Evolution of the Constants of Motion -----**/
    compute_ELCQ_dot(P_p, ECCENTRICITY, INCLINATION,LZ, Q_CONSTANT,energy_0e, angular_momentum_z_0e, q_carter_0e,
                    &ENERGY_dot, &LZ_dot, &Q_CONSTANT_dot);
        
/**----- Compute the Changes in the Constants of Motion and the NEW Constants of Motion -----**/
    compute_new_constants_of_motion(time_elapsed, ENERGY_dot, LZ_dot, Q_CONSTANT_dot, &ENERGY, &LZ, &C_CONSTANT, &Q_CONSTANT);
               
/**----- Compute the NEW Orbital Parameters -----**/
    ELCQ_to_peiz(ENERGY, LZ, C_CONSTANT, Q_CONSTANT,
                &P_p,&ECCENTRICITY, &THETA_inc, &INCLINATION,
                &R_peri, &R_apo, &THETA_min, &Z_minus, 
                &R_3, &R_4, &P3_p, &P4_p, &Z_plus);
        
/**----- Check whether the System has Plunged or not (GR Criterium) -----**/
    check_if_orbit_plunges(ECCENTRICITY,Q_CONSTANT);

/**----- Computing the Keplerian period  -----**/
    compute_Keplerian_period(P_p,ECCENTRICITY,&KEPLER_period,&StepKT);
   
    n_kepler_max = N_periods*StepKT;

/**----- Resetting the Time Elapsed  -----**/
    time_elapsed = 0.0;
  }
/*----- Updating Old Values of the ODE Angles -----*/
    CHI_OLD = CHI_NEW;
    PSI_OLD = PSI_NEW;
    PHI_OLD = PHI_NEW;
     
/***----- Storing the SCO Evolution-----***/
    fprintf(data,"%lu  %4.6e  %4.6e  %4.6e\n",n_global,R_p_o, COSTH_p_o, PHI_NEW);
    
/***----- Checking wheter we have reached separatrix-----***/
      if(R_peri == R_3) { printf("[!] We have hit the separatrix"); exit(0); }
  }
    
    
/*======= EVOLUTION  LOOP  ENDS  HERE! ========*/
    
/*----- This is the end -----*/  
  return 0;
}




/*---------------------------------------------------------------------------------------*/
/*   The equations for psi and chi are also given in equations (A16) and (A3),           */
/*   respectively, of the following reference:                                           */
/*                                                                                       */
/*       S Drasco & S A Hughes, "Rotating black hole orbit functionals in the            */
/*          frequency domain", PHYSICAL REVIEW D 69, 044015 (2004)                       */
/*                                                                                       */
/*            GEODESIC EQUATIONS IN BOYER-LINDQUIST COORDINATES                          */
/*                                                                                       */
/* @ The geodesic equation for the Boyer-Lindquist radial coordinate r is:               */
/*                                                                                       */
/*         4   / dr \ 2        2   2        2          2           2                     */
/*      rho   | ---- |   = [E(r + a ) - a L] - Delta [r + (L - a E)  + Q]                */
/*             \ dT /                                                                    */
/*                                                                                       */
/*   where T denotes proper time and rho and Delta are given by                          */
/*                                                                                       */
/*                        2    2    2    2                                               */
/*                     rho  = r  + a  cos (theta)                                        */
/*                                                                                       */
/*                               2            2                                          */
/*                      Delta = r  - 2 M r + a                                           */
/*                                                                                       */
/*                                                                                       */
/* @ The geodesic equation for the Boyer-Lindquist polar angular coordinate              */
/*   theta is given by:                                                                  */
/*                                                                                       */
/*         4   / d theta \ 2            2         2    2   2            2                */
/*      rho   | --------- |   =  Q - cot (theta) L  - a cos (theta) (1-E )               */
/*             \    dT   /                                                               */
/*                                                                                       */
/*                                                                                       */
/* @ The geodesic equation for the Boyer-Lindquist time coordinate is:                   */
/*                                                                                       */
/*                          2   2 2                                 2   2                */
/*      2  / dt \         (r + a )     2   2                       r + a                 */
/*   rho  | ---- | =  E [---------- - a sin (theta) ] + a L ( 1 - -------  )             */
/*         \ dT /           Delta                                  Delta                 */
/*                                                                                       */
/*                                                                                       */
/* @ The geodesic equation for the Boyer-Lindquist azimuthal angular                     */
/*   coordinate phi is given by:                                                         */
/*                                                2   2             2                    */
/*      2  / d phi \           2                 r + a             a L                   */
/*   rho  | ------- | = L cosec (theta) + a E ( -------- - 1 ) - -------                 */
/*         \  dT   /                             Delta            Delta                  */
/*                                                                                       */
/*                                                                                       */
/*---------------------------------------------------------------------------------------*/
/*                                                                                       */
/*    TRANSFORMATIONS TO AVOID THE TURNING POINTS IN THE GEODESIC EQUATIONS              */
/*                                                                                       */
/* @ Transformation to avoid the turning points in the radial coordinate:                */
/*                                                                                       */
/*                                   p M                                                 */
/*                         r = ----------------                                          */
/*                              1 + e cos(psi)                                           */
/*                                                                                       */
/*   where p is the semi-latus rectum of the orbit and e is the eccentricity             */
/*                                                                                       */
/*                                                                                       */
/* @ Transformations to avoid the turning points in the polar angular coordinate:        */
/*                                                                                       */
/*                                  2                                                    */
/*                           z = cos (theta)                                             */
/*                                                                                       */
/*                   cos(theta) = sqrt(z_) cos(chi)                                      */
/*                                                                                       */
/*   where z_ is the smallest turning point of z (there are two).                        */
/*                                                                                       */
/*****************************************************************************************/












