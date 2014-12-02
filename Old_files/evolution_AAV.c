/*==============================================================================
 
 This routine evolves the geodesic equations for the Kerr Spacetime using 
 osculating evolution conditions (see my notes) together with an action angle 
 variable formulation.
 
 NOTE: Here we are effectively computing the:

 Dissipative SF: Gives the secular change on E, Lz and Q. Mimics the RR effects 
 by its time evolution.
 
 Conservative SF: Shiftes the instantaneous values of the orbital parameters and
 changes the orbital phases
 
 -------------------------------------------------------------------------------
                                            Created by Priscilla on 20/12/2012
 -------------------------------------------------------------------------------
 Last Update: 29.12.12
 =============================================================================*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/types.h>
#include<time.h>

#include "global_quantities_kerr.h"
#include "global_quantities_main.h"
#include "global_quantities_osc.h"
#include "physical_quantities.h"
#include "macros_2.h"

int AAV_BLcoordinates(double, double, double *, double *);
int check_if_orbit_plunges(double, double);
int compute_frequencies(double *, double *, double *, double *, double *,double *, double *, double *, double *);
int compute_Keplerian_period(double, double, double *, long unsigned *);
int compute_new_constants_of_motion(double, double, double, double, double *, double *, double *, double *);
int ELCQ_to_peiz(double, double, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int geodesic_evolution(double *, double *, double *, double, double, double, double, double);
int Mino_time(double, double, double,double, double, double,double, double,double*, double*);
int osculating_evolution(double *, double *, double *,double, double, double,double, double, double, double);
int PN_ELzQdot(double *, double *, double *,double ,double);

int evolution_AAV( )
{
   FILE *data;

  long unsigned i, maxlen;   // Time Step Counter and total evolution
  long unsigned T_tot;       // Total lenght of the evolution
  
  double t_initial, t_in_sec;// evolution time
  double w_th, w_r;          // Action angle variables
  double W_obs_r, W_obs_th;  // Frequencies in Observer time
  double phi_o, phi;         // BL variable phi
  double Msec, dt;           // MBH mass in seconds

  double lambda_r, lambda_theta;

/*-----File to store data-----*/
  char evolution_AAV[100];
  strcpy(evolution_AAV, GV_Dn_run);
  strcat(evolution_AAV, "/Evolution_AAV.txt");
  if(Store_evolution == 'y')
    data = fopen(evolution_AAV,"w");

/*--- Initializing some varibles --*/
  Msec = MASS_in_MSUN*MASS_SUN*(GN/CSPEED3);
  dt = 1.;
  t_in_sec = 0;
  T_tot = ceil(TIME_yrs*(3600*24*365)) ;
  maxlen = T_tot/dt;
 
/*======= EVOLUTION  LOOP  STARTS  HERE! ========*/

/*----- Orbital initial position of the Particle -----*/
  if ( FLAG_equatorial_orbit == 'y' )
  {
    PSI_OLD = PSI_o;  CHI_OLD = 0.5*PI;  PHI_OLD = PHI_o;
  }
  else
  {
    PSI_OLD = PSI_o;  CHI_OLD = CHI_o;  PHI_OLD = PHI_o;
  }
   
  for (i = 0; i <= maxlen; i++)
  {
/*----- Computational times -----*/

    t_in_sec += dt;
      
/*----- Setting the Time, Time Intervals, and Variables for ODE Evolution of the EMRI trajectory -----*/
    t_initial = t_in_sec - DT_MBH*BHMASS_to_s;
	  
    DT_INI = 0.1*DT_MBH*BHMASS_to_s;
      
    CHI_NEW = CHI_OLD;
    PSI_NEW = PSI_OLD;
    PHI_NEW = PHI_OLD;
       
/*---Fundamental frequencies seen by a distant observer---*/
    W_obs_r  = (OMEGA_r/OMEGA_t)*Msec/(2*PI);
    W_obs_th = (OMEGA_th/OMEGA_t)*Msec/(2*PI);
      
/*--- Computing the Action Angle variables wr and wth----*/
    Mino_time(CHI_NEW, CHI_NEW, ENERGY, ECCENTRICITY, P_p, P3_p, P4_p, SPIN, &lambda_r, &lambda_theta);
      
    w_r  = OMEGA_r*lambda_r;
    w_th = OMEGA_th*lambda_theta;

/*---- Analytical expressions for the radial and polar BL coordinates in terms of the AAV ---*/
    AAV_BLcoordinates(w_r, w_th, &R_p_o, &COSTH_p_o);
      
    // printf("\n dt= %4.6e  Lr=%4.6e  Lth=%4.6e  wr=%4.6e  wth=%4.6e   r=%4.6e  costh=%4.6e",dt,lambda_r, lambda_theta,w_r_o,w_th_o, R_p_o, COSTH_p_o);
  
/*---- Evolution along a given geodesic ---*/
    if ( MAG_SF == 0)
        geodesic_evolution(&CHI_NEW, &PSI_NEW, &PHI_NEW, t_initial, t_in_sec, DT_INI, DT_MIN*Msec, ACCURACY);//, &n_ok, &n_bad);
/*---- Evolution through geodesics using Osculating evolution and PN SF effects ---*/
    else 
    {
/*---- Evolution of the constants of motion using PN self-force ---*/
     PN_ELzQdot(&ENERGY_dot, &LZ_dot, &Q_CONSTANT_dot, R_p_o, COSTH_p_o);

/*---- Evolution of the inital phases (positional elements) ---*/
     osculating_evolution(&CHI_NEW, &PSI_NEW, &PHI_NEW,t_initial, t_in_sec, DT_INI, DT_MIN*Msec, ACCURACY,lambda_r,lambda_theta);

/*-----Updating Constants of Motion-----**/
     compute_new_constants_of_motion(dt,ENERGY_dot,LZ_dot,Q_CONSTANT_dot, &ENERGY, &LZ, &C_CONSTANT, &Q_CONSTANT);
          
/*-----Updating Orbital Parameters and Turning points -----*/ 
     ELCQ_to_peiz(ENERGY, LZ, C_CONSTANT, Q_CONSTANT,&P_p, &ECCENTRICITY, &THETA_inc, &INCLINATION, &R_peri, &R_apo, &THETA_min, &Z_minus,&R_3, &R_4, &P3_p, &P4_p, &Z_plus);
        
/**----- Check whether the System has Plunged or not (GR Criterium) -----**/
     check_if_orbit_plunges(ECCENTRICITY,Q_CONSTANT);
  
/*---- Updating fundamental frequencies (associated with Mino time) and Elliptic function arguments ---*/
       compute_frequencies(&OMEGA_r, &OMEGA_th, &OMEGA_t, &Kcom_r, &Kcom_th, &Ecom_r, &Ecom_th, &KR, &KTH);
    }
/*----- Updating Old Values of the ODE Angles -----*/
      CHI_OLD = CHI_NEW;
      PSI_OLD = PSI_NEW;
      PHI_OLD = PHI_NEW;
      
    if (Store_evolution == 'y')
    {
     //fprintf(data,"%4.6e  %4.6e  %4.6e  %4.6e  %4.6e   %4.6e  %4.6e",dt,lambda_r, lambda_theta,w_r,w_th, R_p_o, COSTH_p_o);
     fprintf(data,"%4.6e  %4.6e %4.6e  %4.6e  %4.6e  %4.6e  %4.6e  \n",P_p,ECCENTRICITY, THETA_inc, ENERGY, LZ, Q_CONSTANT, INCLINATION);
    }
      
          
  }

  if (Store_evolution == 'y')
   fclose(data);
    
/*======= EVOLUTION  LOOP  ENDS  HERE! ========*/
    
/*----- This is the end -----*/  
  return 0;
}
