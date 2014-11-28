/*==============================================================================
 
        This routine initializes the variables used in the current run
 
 -------------------------------------------------------------------------------
                                             Created by Priscilla on 20/12/2012
 -------------------------------------------------------------------------------
 Last Update: 14.08.13
 =============================================================================*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "physical_quantities.h"
#include "global_quantities_osc.h"

#include "macros_2.h" 

int checks_errors(int);
int compute_Keplerian_period(double, double, double *, long unsigned, double);
int pez_to_ELCQ(char, double, double, double,double, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int turning_points(double, double, double, double, double, double, double *, double *, double *, double *, double *);
int compute_frequencies(double *, double *, double *, double *, double *,double *, double *, double *, double *);


int initialization(void)
{
  
/*----- Setting Global System Parameters for this Run -----*/
  MASS_in_MSUN = EMRI_Parameter[1];          // Black Hole Mass (in units of the Mass of the Sun)
  SPIN         = EMRI_Parameter[2];          // Black Hole Spin (in units of the Black Hole Mass)
  MASS_SCO     = EMRI_Parameter[3];          // Mass of the Particle orbiting the Black Hole (in units of the Black Hole Mass)
  ECCENTRICITY = EMRI_Parameter[4];          // Initial Eccentricity
  P_p          = EMRI_Parameter[5];          // Initial Semilatus Rectum
  THETA_inc    = EMRI_Parameter[6]*PI/180.0; // Initial Inclination of the Orbit in radians
  THETA_Source = EMRI_Parameter[7];          // EMRI Polar Angle with respect to the Ecliptic Reference Frame.
  PHI_Source   = EMRI_Parameter[8];          // EMRI Azimuthal Angle with respect to the Ecliptic Reference Frame.
  THETA_Kerr   = EMRI_Parameter[9];          // Spin Vector Polar Angle with respect to the Ecliptic Reference Frame.
  PHI_Kerr     = EMRI_Parameter[10];         // Spin Vector Azimuthal Angle with respect to the Ecliptic Reference Frame.
  DL           = EMRI_Parameter[11];	         // Luminosity Distance divided by the SCO mass.
  PSI_o        = EMRI_Parameter[12];         // Boyer-Lindquist Angle psi [associated with r_BL] at t = t_o: psi_o.
  CHI_o        = EMRI_Parameter[13];         // Boyer-Lindquist Angle chi [associated with theta_BL] at t = t_o: chi_o.
  PHI_o        = EMRI_Parameter[14];         // Boyer-Lindquist Angle varphi [identical to varphi_BL] at t = t_o: varphi_o.

/*----- Setting Global computational Parameters for this Run -----*/
  EVAL_arg = 0;                             // Counter to track whether the elliptical integrals have been passed invalid arguments
  DT_MBH = 10.0*DT_INI;                     // ODE Global Time Step
  N_periods = 1.;                           // Number of "periods" between evolving a geodesic

   // printf("NSAMPLING=%4.6e DT_INI=%4.6e DT_MIN=%lu Nfreq=%lu GV_N_max=%lu\n",NSAMPLING, DT_INI,DT_MIN,Nfreq, NMAX);exit(0);
    
/*---- Setting Flags for this Run -----*/
 if ( SPIN < 1.0e-14 )
  {
    FLAG_Schwarzschild = 'y';       // Flag that controls whether the Black Holes is described by the Schwarzschild metric tensor.  VALUES: y (yes), n (no)
	SPIN = 0.0;
  }
  else FLAG_Schwarzschild = 'n';

  if (FLAG_Schwarzschild == 'y')
  {
    Q_CONSTANT= 0.0;
    
    FLAG_equatorial_orbit = 'y';
    THETA_inc = 0.0;
  }


  if ( ECCENTRICITY < 1.0e-14 )
  {
    FLAG_circular_orbit = 'y';      // Flag that controls whether the orbit is circular
	  ECCENTRICITY = 0.0;
  }
  else FLAG_circular_orbit = 'n';
  
  if ( THETA_inc < 1.0e-14 )
  {
    FLAG_equatorial_orbit = 'y';   // Flag that controls whether the orbit is equatorial or not
  	THETA_inc = 0.0;
  }
  else FLAG_equatorial_orbit = 'n';
  
  if ( THETA_inc < 0.0 )           // Flag that controls whether the orbit is prograde or retrograde
    FLAG_orbit_spin = 'r';
  else 
    FLAG_orbit_spin = 'p';           
/**----- Inclination Angle must be in the interval [-PI/2,PI/2] -----**/
    if ( (THETA_inc < -0.5*PI) || (THETA_inc > 0.5*PI) )  checks_errors(3);

//----Orbital Phases initialization
  if(FLAG_RR=='y')
  {
    if(FLAG_equatorial_orbit == 'y' )
    {
      PSI_OLD = PSI_o;
      CHI_OLD = 0.5*PI;
      PHI_OLD = PHI_o;
    }
    else
    {
      PSI_OLD = PSI_o;
      CHI_OLD = CHI_o;
      PHI_OLD = PHI_o;
    }
  }
  else if(FLAG_OS=='y')
  {
    PSI_OLD = PSI_o;
    CHI_OLD = CHI_o;
    PHI_OLD = PHI_o;
    
    TAU = 0.0;   // Proper time
    TIME = 0.0;  // coordinate time
    Qr  = 0.;
    Qth = PI/2;
    ODEcount = 0;     
  }
  
/*----- Global Black Hole Parameters [M^2,S^2,R_+,R_-,sqrt(M^2-a^2)] -----*/
  MASS_ratio = MASS_SCO;                    // mu = m/M
  MASS_ratio2 = SQR(MASS_ratio);            // mu^2
  SPIN2  = SQR(SPIN);                       // (a/M)^2
  R_plus = 1.0 + sqrt(1.0 - SPIN2);          // r+ = M + sqrt(M^2 - a^2)
  R_minus= 1.0 - sqrt(1.0 - SPIN2);          // r- = M - sqrt(M^2 - a^2)
 
/*-----Computing Inital Constants of Motion -----*/
  pez_to_ELCQ(FLAG_circular_orbit, P_p, ECCENTRICITY, THETA_inc,SPIN,
              &INCLINATION, &R_peri, &R_apo, &THETA_min, &Z_minus,
              &ENERGY, &LZ, &C_CONSTANT, &Q_CONSTANT);
    
/*-----Computing initial turning points of the Boyer-Linquist radial coordinate 'r' and 'theta' -----*/
  turning_points(P_p, ECCENTRICITY, ENERGY, LZ,Q_CONSTANT,Z_minus, &R_3, &R_4, &P3_p, &P4_p, &Z_plus);

/*----Computing Initial fundamental frequencies (associated with Mino time) and initial arguments of the Elliptic functions ---*/
  compute_frequencies(&OMEGA_r, &OMEGA_th, &OMEGA_t, &Kcom_r, &Kcom_th, &Ecom_r, &Ecom_th, &KR, &KTH);

  // analytic_r_costh(Qr, Qth, &R_p_o, &COSTH_p_o);
  //Z_minus = COSTH_p_o;  // this horrible think is because Z- is not exactelly costh for chi =0. I found a typo in Fujitas's paper

  // printf("E =%4.6e Lz =%4.6e Q =%4.6e\n", ENERGY, LZ, Q_CONSTANT);
  
/*----- Computing Keplerian period and number of time steps inside it (used here to apply RR) -----*/
  if(FLAG_RR=='y')
   compute_Keplerian_period(P_p,ECCENTRICITY,&KEPLER_period,&StepKT,DT_MBH);
  if ( (FLAG_RR == 'y') && ((StepKT*N_periods) > NMAX) ) checks_errors(23);

    
/*----- This is the end -----*/
  return 0;
}

