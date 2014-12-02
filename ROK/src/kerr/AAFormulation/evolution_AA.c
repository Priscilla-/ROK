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
 Last Update: 30.09.13
 =============================================================================*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/types.h>
#include<time.h>
#include <math.h> 

#include "global_quantities_kerr.h"
#include "global_quantities_main.h"
#include "global_quantities_osc.h"
#include "physical_quantities.h"
#include "macros_2.h" 

int analytic_r_costh(double*,double *,double *);
int check_if_orbit_plunges(double, double);
int compute_frequencies(double *, double *, double *, double *, double *,double *, double *, double *, double *);
int compute_Keplerian_period(double, double, double *, long unsigned *, double);
int evolution_step(double *,double, double, double, double, double );
int geodesic_evolution(double *, double *, double *,double *, double, double, double, double, double);
int waveform( double *, double* , double * ,double , double);

double *Realvector(const long nl, const long nh);
void free_Realvector(double *v, const long nl, const long nh);

int evolution_AA(FILE *data, FILE *dataf, FILE *dataflux,FILE *dataGW, double *x)
{
   double step_insp;
   double dt_sec, t_evol;                                    // vector of variables
   double r, costh;
   double h_p, h_x;
   int i;

//----Time interval for data saving
  dt_sec = DTsec;//STORE_rate;
  step_insp = (double) (TIME_yrs*SECSYR/dt_sec);
  if (step_insp < 3) step_insp = 3;

  //printf("%4.6e\n",step_insp); exit(0);
  
/*======= EVOLUTION  LOOP  STARTS  HERE! ========*/
  t_evol = 0.;
   
  for (i = 0; i <= step_insp; i++)
  {
    if(MAG_SF != 0)
    {
      x[1] = PHI_OLD;
      x[2] = TAU;
      x[3] = P_p;
      x[4] = ECCENTRICITY;
      x[5] = THETA_inc;
      x[6] = Qr;
      x[7] = Qth;
      x[8] = Lambda;
     }
    else
    {
      x[1] = Qr;
      x[2] = Qth;
      x[3] = PHI_OLD;
      x[4] = TIME;
    }
    
/*-- Computing analytical expressions for r(w_r) and cos(th)[w_th] and the radial and polar phases psi & chi --*/
    analytic_r_costh(x, &r, &costh);

    // printf("%4.12e  %4.12e  %4.12e %4.12e  %4.12e  %4.12e %4.12e\n",t_evol, r, costh, COS_PSI, COS_CHI, Qr, Qth);
    // printf("Tau=%4.6e, P=%4.6e,   Ecc=%4.6e, theta_inc=%4.6e\n", TAU, P_p,ECCENTRICITY, THETA_inc);
    // printf("E=%4.6e Lz =%4.6e, Q = %4.6e\n\n",ENERGY,LZ, Q_CONSTANT);
  
/*--- ODEs evolution--*/
    evolution_step(x, t_evol, t_evol+dt_sec,DT_MIN/BHMASS_to_s,dt_sec, ACCURACY/BHMASS_to_s);
    
/*----- Compute Gravitational Waveforms and Response Functions -----*/
    waveform(&h_p, &h_x, x, r, costh);
    
    //printf("h_p=%4.6e h_x =%4.6e\n",h_p,h_x);
    
    if(MAG_SF != 0)
    {
      fprintf(dataflux,"%4.6e  %4.12e  %4.12e  %4.14e\n",t_evol, ENERGY_dot,LZ_dot,Q_CONSTANT_dot);
      fflush(dataflux);
    }
      
    fprintf(dataf   ,"%4.12e  %4.12e  %4.12e %4.12e  %4.12e  %4.12e %4.12e\n",t_evol, r, costh, COS_PSI, COS_CHI, Qr, Qth);
    fflush(dataf);
      
    fprintf(data    ,"%4.6e  %4.6e  %4.6e  %4.6e  %4.6e  %4.6e  %4.6e  %4.6e\n",t_evol, P_p,ECCENTRICITY, THETA_inc, ENERGY, LZ, Q_CONSTANT,INCLINATION);
    fflush(data);
      
    fprintf(dataGW  ,"%4.6e  %4.6e  %4.6e  \n",t_evol, h_p, h_x);
    fflush(dataGW);
  
/**----- Check whether the system has plunged -----**/
    check_if_orbit_plunges(ECCENTRICITY,Q_CONSTANT);
    
/*----- Updating Old Values -----*/    
    if(MAG_SF != 0)
    {
     PHI_OLD = x[1];
     TAU = x[2];
     P_p = x[3];
     ECCENTRICITY = x[4];
     THETA_inc  = x[5];
     Qr  = x[6];
     Qth = x[7];
     Lambda = x[8];
    }
    else
    {
     Qr  = x[1];
     Qth = x[2];
     PHI_OLD= x[3];
     TIME = x[4];
    }
  
   t_evol += dt_sec;
  }
/*======= EVOLUTION  LOOP  ENDS  HERE! ========*/
    
/*----- This is the end -----*/  
  return 0;
}
