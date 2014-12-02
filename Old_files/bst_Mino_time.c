//
//  bst_Mino_time.c
//  X_ROK_V1.7
//
//  Created by Priscilla on 24/07/2013.
//
//


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"

#include "macros_1.h"

#define K_MAXIMUM 8      //  Maximum row number in the extrapolation
#define KSEC (K_MAXIMUM+2)   //
#define SAFETY_1 0.25      //  Safety factors
#define SAFETY_2 0.7       //
#define TINY 1.0e-30       //  Prevents division by zero
#define RED_MIN 0.7      //  Minimum factor for step size reduction
#define RED_MAX 1.0e-5     //  Maximum factor for step size reduction
#define SCALE_MAX 0.1      //  1/SCALE_MAX = Maximum factor for step increasing

int rhs_Mino_time(double *, double , double );
//Defined below
int midpoint_method(long unsigned, double, double, double *);
int poly_extrapol(long unsigned, double, double, double, double, double *, double *, double *, double *, double *, double *);

double tableau[3][KSEC],zz[KSEC];


bst_Mino_time(&lambda,&d_initial,dstep,accuracy,lamda_scale,&dstep_next);


int bst_Mino_time(double *lambda, double *d_initial, double dt_try, double tolerance, double lambda_scale, double *dt_next)
{
  double dt_used[3];
  double lmbd_save[3];         // Variables to save Chi_p, Psi_p, and Phi_p
  double lmbd_estimated[3];    // Values of chi, psi, and phi from the midpoint method
  double lmbd_extrapolated;    // Values of chi, psi, and phi after extrapolation
  double lmbd_err;             // Estimated errors from extrapolation
  double tolerance_1;          // Tolerances
  
  double dstep[3];             // Step
  double z_ind[3];             // Independent variable for polynomial extrapolation
  double error_max;            // Normalized error estimate
  double factor,reduction_factor,scale;
  double work,work_min;
  double error[KSEC];         // Normalized error estimate
    
  static double time_new[3];
  static double old_tolerance=-1.0;
  static double a[KSEC];           // Work coefficients
  static double alpha[KSEC][KSEC];     // alpha(k,q)
  static long unsigned k_max,k_optimal;  // Optimal row number for convergence
  static long unsigned first=1;
    
  long unsigned k,kk,iq,km;
  long unsigned reduction;
  long unsigned exit_flag=0;
  
  lambda_estimated[1] = 0.0;
  lambda_estimated[2] = 0.0;
  lambda_estimated[3] = 0.0;
  
  if (tolerance != old_tolerance)
  {
    *dt_next[1] = time_new[1] = -1.0e29;
    *dt_next[2] = time_new[2] = -1.0e29;
    *dt_next[3] = time_new[3] = -1.0e29;
      
    tolerance_1 = SAFETY_1*tolerance;
    a[1] = 3;
    for (k=1;k<=K_MAXIMUM;k++)
    {
      a[k+1] = a[k] + 2*(k+1);
    }
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
  
  dstep = dt_try;
  
  lmbd_save = *lambda;

  
  if (( *d_initial[1] != time_new[1] || dstep[1] != (*dt_next[1]) )&&( *d_initial[2] != time_new[2] || dstep[2] != (*dt_next[2]) )&&( *d_initial[3] != time_new[3] || dstep[3] != (*dt_next[3]) ))
  {
    first = 1;
    k_optimal = k_max;
  }
  
  reduction = 0;
  
  for(;;)
  {
    for (k=1;k<=k_max;k++)
    {
      time_new[1] = (*d_initial[1]) + dstep[1];
      if (time_new[1] == (*d_initial[1])){ printf(" Lambda(time) | Step size underflow!!!\n");exit(0);}
      
      time_new[2] = (*d_initial[2]) + dstep[2];
        if (time_new[2] == (*d_initial[2])){ printf(" Lambda(psi) | Step size underflow!!!\n");exit(0);}
  
      time_new[3] = (*d_initial[3]) + dstep[3];
        if (time_new[3] == (*d_initial[3])){ printf(" Lambda(psi) | Step size underflow!!!\n");exit(0);}
        
      midpoint_method(k,dstep,lmbd_save,&lmbd_estimated);
      
      z_ind[1] = dstep[1]*dstep[1]/k/k/4.0;
      z_ind[2] = dstep[2]*dstep[2]/k/k/4.0;
      z_ind[3] = dstep[3]*dstep[3]/k/k/4.0;
        
      poly_extrapol(k,z_ind,lmbd_estimated,&lmbda_extrapolated,&lmbd_err);
      
      *lambda = lmbd_extrapolated;
            
      if (k != 1)
      {
        error_max = TINY;
        error_max = DMAX(error_max,fabs(lmbd_err[1]/lambda_scale[1]));
        error_max = DMAX(error_max,fabs(lmbd_err[2]/lambda_scale[2]));
        error_max = DMAX(error_max,fabs(lmbd_err[3]/lambda_scale[3]));
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
    dstep[1] *= reduction_factor;
    dstep[2] *= reduction_factor;
    dstep[3] *= reduction_factor;
    reduction = 1;
  }
  
  *d_initial = time_new;
   dt_used = dstep;
    
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
  
  *dt_next = dstep/scale;
  
  if (k_optimal>=k && k_optimal != k_max && !reduction)
  {
    factor = DMAX(scale/alpha[k_optimal-1][k_optimal],SCALE_MAX);
    if (a[k_optimal+1]*factor <= work_min)
    {
      *dt_next = dstep/factor;
      k_optimal++;
    }
  }
  
  return 0;
}


int midpoint_method(long unsigned n_half_steps, double dstep_total, double lambda, double *lmbd_predicted)
{
  long unsigned n;
  double z_minus_lmbd[3];
  double dz_minus_lmbd[3];
  double z_plus_lmbd[3];
  double dz_plus_chi,dz_plus_psi,dz_plus_phi;
  double dz_store_chi,dz_store_psi,dz_store_phi;
  double derivative_chi,derivative_psi,derivative_phi;
  double dstep[3],dt_2;
  
  
  /*----- Step Size for this iteration -----*/
  dstep = dstep_total/(2*n_half_steps);
  
  
  /*----- ZEROth Stage -----*/
  z_minus_lmbd = lambda;
 
  
  dz_minus_lmbd[1] = 0.0;
  dz_minus_lmbd[2] = 0.0;
  dz_minus_lmbd[3] = 0.0;
  
  /*----- FIRST Stage -----*/
   rhs_Mino_time(&derivative_lmbd,psi, chi);
    
  dz_plus_lmbd = dstep*derivative_lmbd;
  
  z_plus_lmbd = z_minus_lmbd + dz_plus_lmbd;
 
  
  
  /*----- GENERAL Step -----*/
  if(FLAG_RR=='y')
    geodesic_RHS(&derivative_chi,&derivative_psi,&derivative_phi,z_plus_chi,z_plus_psi);
  else
    osculating_RHS(&derivative_chi,&derivative_psi,&derivative_phi,z_plus_chi,z_plus_psi);
  
  dt_2 = 2.0*dt;
  
  for (n=2;n<=2*n_half_steps;n++)
  {
    dz_store_chi = dz_minus_chi + dt_2*derivative_chi;
    dz_store_psi = dz_minus_psi + dt_2*derivative_psi;
    dz_store_phi = dz_minus_phi + dt_2*derivative_phi;
    
    dz_minus_chi = dz_plus_chi;
    dz_minus_psi = dz_plus_psi;
    dz_minus_phi = dz_plus_phi;
    
    dz_plus_chi = dz_store_chi;
    dz_plus_psi = dz_store_psi;
    dz_plus_phi = dz_store_phi;
    
    z_plus_chi = z_minus_chi + dz_plus_chi;
    z_plus_psi = z_minus_psi + dz_plus_psi;
    /*   z_plus_phi = z_minus_phi + dz_plus_phi; */
    
    
    if(FLAG_RR=='y')
      geodesic_RHS(&derivative_chi,&derivative_psi,&derivative_phi,z_plus_chi,z_plus_psi);
    else
      osculating_RHS(&derivative_chi,&derivative_psi,&derivative_phi,z_plus_chi,z_plus_psi);
    
  }
  
  dz_store_chi = 0.5*(dz_minus_chi+dz_plus_chi+dt*derivative_chi);
  dz_store_psi = 0.5*(dz_minus_psi+dz_plus_psi+dt*derivative_psi);
  dz_store_phi = 0.5*(dz_minus_phi+dz_plus_phi+dt*derivative_phi);
  
  *chi_predicted = z_minus_chi + dz_store_chi;
  *psi_predicted = z_minus_psi + dz_store_psi;
  *phi_predicted = z_minus_phi + dz_store_phi;
  
  return 0;
}


int poly_extrapol(long unsigned call_number,double zi,double chi_mid, double psi_mid, double phi_mid, double *chi_extrapol, double *psi_extrapol, double *phi_extrapol, double *dchi, double *dpsi, double *dphi)
{
  double c_chi,c_psi,c_phi;
  double delta_chi,delta_psi,delta_phi;
  double q_chi,q_psi,q_phi;
  double delta,factor_1,factor_2;
  long unsigned j;
  
  zz[call_number] = zi;
  
  *chi_extrapol = chi_mid;
  *psi_extrapol = psi_mid;
  *phi_extrapol = phi_mid;
  
  *dchi = chi_mid;
  *dpsi = psi_mid;
  *dphi = phi_mid;
  
  if (call_number==1)
  {
    tableau[1][1] = chi_mid;
    tableau[2][1] = psi_mid;
    tableau[3][1] = phi_mid;
  }
  else
  {
    c_chi = chi_mid;
    c_psi = psi_mid;
    c_phi = phi_mid;
    
    for (j=1;j<call_number;j++)
    {
      delta = 1.0/(zz[call_number-j]-zi);
      factor_1 = zi*delta;
      factor_2 = zz[call_number-j]*delta;
      
      q_chi = tableau[1][j];
      q_psi = tableau[2][j];
      q_phi = tableau[3][j];
      
      tableau[1][j] = *dchi;
      tableau[2][j] = *dpsi;
      tableau[3][j] = *dphi;
      
      delta_chi = c_chi - q_chi;
      delta_psi = c_psi - q_psi;
      delta_phi = c_phi - q_phi;
      
      *dchi = factor_1*delta_chi;
      *dpsi = factor_1*delta_psi;
      *dphi = factor_1*delta_phi;
      
      c_chi = factor_2*delta_chi;
      c_psi = factor_2*delta_psi;
      c_phi = factor_2*delta_phi;
      
      *chi_extrapol += *dchi;
      *psi_extrapol += *dpsi;
      *phi_extrapol += *dphi;
    }
    
    tableau[1][call_number] = *dchi;
    tableau[2][call_number] = *dpsi;
    tableau[2][call_number] = *dphi;
  }
  
  
  return 0;
  
}
