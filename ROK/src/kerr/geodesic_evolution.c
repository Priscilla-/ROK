/*==============================================================================
 
  This routine evolves the geodesic equation in the KERR spacetime to get the
  new particle's position and velocity.  These equations are ODEs. The numerical
  scheme used for the integration is:
 
                  Bulirsch-Stoer (See numerical recipes in C).
 
--------------------------------------------------------------------------------
 Last Update : 08.07.13                                                               
==============================================================================*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"

#include "physical_quantities.h"
#include "macros_2.h"

#define TINY 1.0e-30
#define MAX_NUM_STEPS 10000

int geodesic_RHS(double *, double *, double *,double *, double, double);
int bulirsch_stoer(double *, double *,double *, double *, double *, double, double, double, double, double, double, double *, double *);

int geodesic_evolution(double *psi, double *chi, double *phi,double *lambda, double t_initial, double t_final, double dt_initial, double dt_min, double accuracy)
{
    
    
  double chi_save,psi_save,phi_save,lambda_save;               // Variables to save chi, psi, and phi
  double dchi_dt,dpsi_dt,dphi_dt,dlambda_dt;                  // Derivatives of chi, psi, and phi
  double chi_scale,psi_scale,phi_scale, lambda_scale;            // Scaling of the errors
  double time,dt;                                  // Time and Time-step size                 
  double dt_used,dt_next;                          // Time step used by the BS routine and the one computed for the next step      
  long unsigned n;
  
  time = t_initial;
  dt = SIGN(dt_initial,t_final-t_initial); // = DT_MBH
 

  chi_save = *chi; 
  psi_save = *psi;
  phi_save = *phi;
  lambda_save = *lambda;
    
 
  for (n=1;n<=MAX_NUM_STEPS;n++)
  {
    geodesic_RHS(&dchi_dt,&dpsi_dt,&dphi_dt,&dlambda_dt,chi_save,psi_save);    /* Computing derivatives (Evaluating the RHSs) */
      
    chi_scale = fabs(chi_save) + fabs(dchi_dt*dt) + TINY;
    psi_scale = fabs(psi_save) + fabs(dpsi_dt*dt) + TINY;
    phi_scale = fabs(phi_save) + fabs(dphi_dt*dt) + TINY;
    lambda_scale = fabs(lambda_save) + fabs(dlambda_dt*dt) + TINY;
      
    //  printf(" psi_scale =%4.6e psi_save =%4.6e, dpsi_dt =%4.6e, dt = %4.6e \n", psi_scale,psi_save,dpsi_dt, dt);
    //  printf(" chi_scale =%4.6e chi_save =%4.6e, dchi_dt =%4.6e, dt = %4.6e \n", chi_scale,chi_save,dchi_dt, dt);
    //  printf(" phi_scale =%4.6e phi_save =%4.6e, dphi_dt =%4.6e, dt = %4.6e \n", phi_scale,phi_save,dphi_dt, dt);
      
      
    if ( (time+dt-t_final)*(time+dt-t_initial) > 0.0 )                /* Decreases the time step if it is going to overshoot */
    {
      dt = t_final - time;
    }

    bulirsch_stoer(&chi_save,&psi_save,&phi_save,&lambda_save,&time,dt,accuracy,chi_scale,psi_scale,phi_scale,lambda_scale,&dt_used,&dt_next);

    if ( (time-t_final)*(t_final-t_initial) >= 0)
    {
      *chi = chi_save;
      *psi = psi_save;
      *phi = phi_save;
      *lambda = lambda_save;
	  
      return 0; 
    }

    if (fabs(dt_next) <= dt_min) printf("  | Step size too small for computing derivatives!!\n");

    dt = dt_next;  
  }

  printf("  | Too many steps!!!\n");exit(0);

  
/*----- This is the end -----*/  
  return 0;

}
