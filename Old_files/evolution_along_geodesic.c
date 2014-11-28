/*==============================================================================
 
  This routine evolves the geodesic equation in the KERR spacetime to get the
  new particle's position and velocity.  These equations are ODEs. The numerical
  scheme used for the integration is:
 
                  Bulirsch-Stoer (See numerical recipes in C).
 
--------------------------------------------------------------------------------
 Last Update : 10.04.13                                                               
==============================================================================*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>


#include "macros_2.h"

#define TINY 1.0e-30
#define MAX_NUM_STEPS 10000

int RHs_r_theta_phi(double*, double *, double *, double, double);
int bulirsch_stoer_step_AAF(double *, double *, double *, double *, double, double, double, double, double, double *, double *);

int evolution_along_geodesic(double *r_p, double *cosTh_p, double *phi_p, double t_initial, double t_final,double dt_initial, double dt_min, double accuracy, long unsigned *n_ok, long unsigned *n_bad)
{
    
  FILE *data;
    
  double r_save,cosTh_save,phi_save;               // Variables to save chi, psi, and phi
  double dr_dl,dcosTh_dl,dphi_dl;                  // Derivatives of r(wr), costh(wth), and phi
  double r_scale,cosTh_scale,phi_scale;              // Scaling of the errors
  double time,dt, dl;                              // Time and Time-step size
  double dt_used,dt_next;                          // Time step used by the BS routine and the one computed for the next step      
  long unsigned n;
    
  //  double dr_dl,dcosTh_dl,dphi_dl;                  // Derivatives of r(wr), costh(wth), and phi
    
  time = t_initial;
  dt = SIGN(dt_initial,t_final-t_initial);
    
  dl = dt;  // check time meaning here!!!!!!!
 
  *n_ok = 0;
  *n_bad = 0;

  r_save = *r_p;
  cosTh_save = *cosTh_p;
  phi_save = *phi_p;
    

  for (n=1;n<=MAX_NUM_STEPS;n++)
  {
     // printf("n %lu \n",n);
     /* Computing derivatives (Evaluating the RHSs) */
    RHs_r_theta_phi(&dr_dl, &dcosTh_dl, &dphi_dl, r_save, cosTh_save);

    r_scale     = fabs(r_save)     + fabs(dr_dl*dt)     + TINY;            // chi_f = Chi_o + dchi
	cosTh_scale = fabs(cosTh_save) + fabs(dcosTh_dl*dt) + TINY;
    phi_scale   = fabs(phi_save)   + fabs(dphi_dl*dt)   + TINY;

    if ( (time+dt-t_final)*(time+dt-t_initial) > 0.0 )                /* Decreases the time step if it is going to overshoot */
      dt = t_final - time;
    

    bulirsch_stoer_step_AAF(&r_save,&cosTh_save,&phi_save,&time,dt,accuracy,r_scale,cosTh_scale,phi_scale,&dt_used,&dt_next);

    if (dt_used == dt)
      ++(*n_ok);
    else
      ++(*n_bad);

    if ( (time-t_final)*(t_final-t_initial) >= 0)
    {
      *r_p     = r_save;
	  *cosTh_p = cosTh_save;
      *phi_p   = phi_save;
	  
      return 0; 
    }

    if (fabs(dt_next) <= dt_min) printf("  | Step size too small for computing derivatives!!\n");

    dt = dt_next;  
  }

  printf("  | Too many steps!!!\n");
 
    
    
/*----- This is the end -----*/  
  return 0;

}
