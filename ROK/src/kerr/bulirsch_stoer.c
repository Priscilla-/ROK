/*===================================================================
  This routine carries out a Bulirsch-Stoer step with monitoring of
  local truncation error to ensure accuracy and adjust stepsize.
 
  This has been adapted from Numerical Recipes in C.
---------------------------------------------------------------------
 Last Update : 29.07.2013
====================================================================*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"

#include "macros_1.h"

#define K_MAXIMUM 8         // Maximum row number in the extrapolation
#define KSEC (K_MAXIMUM+2)  //
#define SAFETY_1 0.25       // Safety factors
#define SAFETY_2 0.7        //
#define TINY 1.0e-30        // Prevents division by zero
#define RED_MIN 0.7         // Minimum factor for step size reduction
#define RED_MAX 1.0e-5      // Maximum factor for step size reduction
#define SCALE_MAX 0.1       // 1/SCALE_MAX = Maximum factor for step increasing

int midpoint_method(long unsigned, double, double, double, double, double, double *,double *, double *, double *);
int geodesic_RHS(double *, double *, double *, double *, double, double);
int poly_extrapol(long unsigned, double, double, double, double, double, double *,double *,double *, double *, double *, double *, double *, double *);

double tableau[4][KSEC],zz[KSEC];

int bulirsch_stoer(double *chi, double *psi, double *phi, double *lambda, double *time, double dt_try, double tolerance, double chi_scale, double psi_scale, double phi_scale, double lambda_scale, double *dt_used, double *dt_next)
{
   double chi_save,psi_save,phi_save, lambda_save;             // Variables to save Chi_p, Psi_p, and Phi_p
   double chi_estimated,psi_estimated,phi_estimated, lambda_estimated;   // Values of chi, psi, and phi from the midpoint method
   double chi_extrapolated,psi_extrapolated,phi_extrapolated, lambda_extrapolated; // Values of chi, psi, and phi after extrapolation
   double chi_err,psi_err,phi_err,lambda_err; // Estimated errors from extrapolation
   double tolerance_1; // Tolerances
   double dt;          // Step
   double z_ind;       // Independent variable for polynomial extrapolation
   double error_max;   // Normalized error estimate
   double factor,reduction_factor,scale;
   double work,work_min;
   double error[KSEC]; // Normalized error estimate
  
   static double a[KSEC];           // Work coefficients
   static double alpha[KSEC][KSEC]; // alpha(k,q)
   static double old_tolerance=-1.0;
  
   static long unsigned k_max,k_optimal; // Optimal row number for convergence
   static long unsigned first=1;
   static double time_new;
 
   long unsigned exit_flag=0;
   long unsigned k,kk,iq,km;
   long unsigned reduction;
  
   chi_estimated = 0.0;
   psi_estimated = 0.0;
   phi_estimated = 0.0;
   lambda_estimated = 0.0;
   
   if (tolerance != old_tolerance)
   {
     *dt_next = time_new = -1.0e29;
     tolerance_1 = SAFETY_1*tolerance;
     a[1] = 3;
     
     for (k=1;k<=K_MAXIMUM;k++)
      a[k+1] = a[k] + 2*(k+1);
     
     for (iq=2;iq<=K_MAXIMUM;iq++)
     {
      for (k=1;k<iq;k++)
      {
      alpha[k][iq] = pow(tolerance_1,(a[k+1] - a[iq+1])/((a[iq+1] - a[1] + 1.0)*(2*k+1)));
      }
     }
     old_tolerance = tolerance;
     for (k_optimal=2;k_optimal<K_MAXIMUM;k_optimal++)
     {
      if (a[k_optimal+1]>a[k_optimal]*alpha[k_optimal-1][k_optimal])
      break;
     }
     k_max = k_optimal;
   }
   
   dt = dt_try;
   
   chi_save = *chi;
   psi_save = *psi;
   phi_save = *phi;
   lambda_save = *lambda;
   
   if ( *time != time_new || dt != (*dt_next) )
   {
     first = 1;
     k_optimal = k_max;
   }
   
   reduction = 0;
   
   for(;;)
   {
    for (k=1;k<=k_max;k++)
    {
      time_new = (*time) + dt;
      if (time_new == (*time)){ printf("  | Step size underflow!!!\n");exit(0);}
      
      midpoint_method(k,dt,chi_save,psi_save,phi_save, lambda_save,&chi_estimated,&psi_estimated,&phi_estimated, &lambda_estimated);
      
      z_ind = dt*dt/k/k/4.0;
      
      poly_extrapol(k,z_ind,chi_estimated,psi_estimated,phi_estimated,lambda_estimated,&chi_extrapolated,&psi_extrapolated,&phi_extrapolated,&lambda_extrapolated,&chi_err,&psi_err,&phi_err,&lambda_err);
      
      *chi = chi_extrapolated;
      *psi = psi_extrapolated;
      *phi = phi_extrapolated;
      *lambda = lambda_extrapolated;
      
      if (k != 1)
      {
      error_max = TINY;
      error_max = DMAX(error_max,fabs(chi_err/chi_scale));
      error_max = DMAX(error_max,fabs(psi_err/psi_scale));
      error_max = DMAX(error_max,fabs(phi_err/phi_scale));
      error_max = DMAX(error_max,fabs(lambda_err/lambda_scale));
      error_max /= tolerance;
      km  = k - 1;
      error[km] = pow(error_max/SAFETY_1,1.0/(2*km+1));
      }
      if ( k != 1 && (k>=(k_optimal-1) || first) )
      {
      if (error_max<1.0)
      {
       exit_flag = 1;
       break;
      }
      
      if ( k == k_max || k == (k_optimal + 1) )
      {
       reduction_factor = SAFETY_2/error[km];
       break;
      }
      else if (k == k_optimal && alpha[k_optimal-1][k_optimal]<error[km])
      {
       reduction_factor = 1.0/error[km];
       break;
      }
      else if (k_optimal == k_max && alpha[km][k_max-1]<error[km])
      {
       reduction_factor = alpha[km][k_max-1]*SAFETY_2/error[km];
       break;
      }
      else if (alpha[km][k_optimal]<error[km])
      {
       reduction_factor = alpha[km][k_optimal-1]/error[km];
       break;
      }
     }
    }
     
     if (exit_flag) break;
     reduction_factor = DMIN(reduction_factor,RED_MIN);
     reduction_factor = DMAX(reduction_factor,RED_MAX);
     dt *= reduction_factor;
     reduction = 1;
   }
   
   *time = time_new;
   *dt_used = dt;
   first = 0;
   work_min = 1.0e35;
   for (kk=1;kk<=km;kk++)
   {
     factor = DMAX(error[kk],SCALE_MAX);
     work = factor*a[kk+1];
     if (work<work_min)
     {
      scale = factor;
      work_min = work;
      k_optimal = kk + 1;
     }
   }
   
   *dt_next = dt/scale;
   
   if (k_optimal>=k && k_optimal != k_max && !reduction)
   {
     factor = DMAX(scale/alpha[k_optimal-1][k_optimal],SCALE_MAX);
     if (a[k_optimal+1]*factor <= work_min)
     {
      *dt_next = dt/factor;
      k_optimal++;
     }
   }
   
   return 0;
}


int midpoint_method(long unsigned n_half_steps, double dt_total, double chi, double psi, double phi,double lambda, double *chi_predicted, double *psi_predicted, double *phi_predicted, double *lambda_predicted)
{
   long unsigned n;
   double z_minus_chi,z_minus_psi,z_minus_phi, z_minus_lambda;
   double dz_minus_chi,dz_minus_psi,dz_minus_phi,dz_minus_lambda;
   double z_plus_chi,z_plus_psi,z_plus_phi,z_plus_lambda;
   double dz_plus_chi,dz_plus_psi,dz_plus_phi, dz_plus_lambda;
   double dz_store_chi,dz_store_psi,dz_store_phi,dz_store_lambda;
   double derivative_chi,derivative_psi,derivative_phi, derivative_lambda;
   double dt,dt_2;
   
   
/*----- Step Size for this iteration -----*/
   dt = dt_total/(2*n_half_steps);
   
/*----- ZEROth Stage -----*/
   z_minus_chi = chi;
   z_minus_psi = psi;
   z_minus_phi = phi;
   z_minus_lambda=lambda;
   
   dz_minus_chi = 0.0;
   dz_minus_psi = 0.0;
   dz_minus_phi = 0.0;
   dz_minus_lambda=0.0;

/*----- FIRST Stage -----*/
    geodesic_RHS(&derivative_chi,&derivative_psi,&derivative_phi,&derivative_lambda,z_minus_chi,z_minus_psi);
    
   dz_plus_chi = dt*derivative_chi;
   dz_plus_psi = dt*derivative_psi;
   dz_plus_phi = dt*derivative_phi;
   dz_plus_lambda=dt*derivative_lambda;
  
   z_plus_chi = z_minus_chi + dz_plus_chi;
   z_plus_psi = z_minus_psi + dz_plus_psi;
   z_plus_phi = z_minus_phi + dz_plus_phi;
   z_plus_lambda = z_minus_lambda + dz_plus_lambda;

/*----- GENERAL Step -----*/
    geodesic_RHS(&derivative_chi,&derivative_psi,&derivative_phi,&derivative_lambda,z_plus_chi,z_plus_psi);
    
   dt_2 = 2.0*dt;
   
   for (n=2;n<=2*n_half_steps;n++)
   {
     dz_store_chi = dz_minus_chi + dt_2*derivative_chi;
     dz_store_psi = dz_minus_psi + dt_2*derivative_psi;
     dz_store_phi = dz_minus_phi + dt_2*derivative_phi;
     dz_store_lambda=dz_minus_lambda + dt_2*derivative_lambda;
     
     
     dz_minus_chi = dz_plus_chi;
     dz_minus_psi = dz_plus_psi;
     dz_minus_phi = dz_plus_phi;
     dz_minus_lambda = dz_plus_lambda;
     
     dz_plus_chi = dz_store_chi;
     dz_plus_psi = dz_store_psi;
     dz_plus_phi = dz_store_phi;
     dz_plus_lambda =dz_store_lambda;
     
     z_plus_chi = z_minus_chi + dz_plus_chi;
     z_plus_psi = z_minus_psi + dz_plus_psi;
     z_plus_phi = z_minus_phi + dz_plus_phi;
     z_plus_lambda = z_minus_lambda+ dz_plus_lambda;
     
     geodesic_RHS(&derivative_chi,&derivative_psi,&derivative_phi,&derivative_lambda,z_plus_chi,z_plus_psi);
  
   }
   
   dz_store_chi = 0.5*(dz_minus_chi+dz_plus_chi+dt*derivative_chi);
   dz_store_psi = 0.5*(dz_minus_psi+dz_plus_psi+dt*derivative_psi);
   dz_store_phi = 0.5*(dz_minus_phi+dz_plus_phi+dt*derivative_phi);
   dz_store_lambda = 0.5*(dz_minus_lambda + dz_plus_lambda + dt*derivative_lambda);
  
   *chi_predicted = z_minus_chi + dz_store_chi;
   *psi_predicted = z_minus_psi + dz_store_psi;
   *phi_predicted = z_minus_phi + dz_store_phi;
   *lambda_predicted = z_minus_lambda + dz_store_lambda;
  
  
   return 0;
}


int poly_extrapol(long unsigned call_number,double zi,double chi_mid, double psi_mid, double phi_mid, double lambda_mid, double *chi_extrapol, double *psi_extrapol, double *phi_extrapol, double *lambda_extrapol, double *dchi, double *dpsi, double *dphi, double *dlambda)
{
   double c_chi,c_psi,c_phi, c_lambda;
   double delta_chi,delta_psi,delta_phi, delta_lambda;
   double q_chi,q_psi,q_phi, q_lambda;
   double delta,factor_1,factor_2;
   long unsigned j;
   
   zz[call_number] = zi;
   
   *chi_extrapol = chi_mid;
   *psi_extrapol = psi_mid;
   *phi_extrapol = phi_mid;
   *lambda_extrapol = phi_mid;
   
   *dchi = chi_mid;
   *dpsi = psi_mid;
   *dphi = phi_mid;
   *dlambda = lambda_mid;
   
   if (call_number==1)
   {
     tableau[1][1] = chi_mid;
     tableau[2][1] = psi_mid;
     tableau[3][1] = phi_mid;
     tableau[4][1] = lambda_mid;
   }
   else
   {
     c_chi = chi_mid;
     c_psi = psi_mid;
     c_phi = phi_mid;
     c_lambda = lambda_mid;
     
     for (j=1;j<call_number;j++)
     {
      delta = 1.0/(zz[call_number-j]-zi);
      factor_1 = zi*delta;
      factor_2 = zz[call_number-j]*delta;
      
      q_chi = tableau[1][j];
      q_psi = tableau[2][j];
      q_phi = tableau[3][j];
      q_lambda = tableau[4][j];
      
      tableau[1][j] = *dchi;
      tableau[2][j] = *dpsi;
      tableau[3][j] = *dphi;
      tableau[4][j] = *dlambda;
       
      delta_chi = c_chi - q_chi;
      delta_psi = c_psi - q_psi;
      delta_phi = c_phi - q_phi;
      delta_lambda = c_lambda-q_lambda;
       
      *dchi = factor_1*delta_chi;
      *dpsi = factor_1*delta_psi;
      *dphi = factor_1*delta_phi;
      *dlambda = factor_1*delta_lambda;
       
      c_chi = factor_2*delta_chi;
      c_psi = factor_2*delta_psi;
      c_phi = factor_2*delta_phi;
      c_lambda=factor_2*delta_lambda;
       
      *chi_extrapol += *dchi;
      *psi_extrapol += *dpsi;
      *phi_extrapol += *dphi;
      *lambda_extrapol +=*dlambda;
     }
     
     tableau[1][call_number] = *dchi;
     tableau[2][call_number] = *dpsi;
     tableau[3][call_number] = *dphi;
     tableau[4][call_number] = *dlambda;
   }
  
   return 0;
}
