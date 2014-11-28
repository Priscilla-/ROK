/*==============================================================================
 This routine computes the analytical expressions for the Mino time lambda,
 parameterized in different ways: Lamda(r) and Lamba(theta), being r and theta 
 BL coordinates
 
 See Fujita and Hikida Class. Quantum Grav. 26 (2009) 135002
 -------------------------------------------------------------------------------
                                            Created by Priscilla on 24/07/2013
 -------------------------------------------------------------------------------
 Last Update: 28/07/2013                                        Modified by PCM
 =============================================================================*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "physical_quantities.h"
#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"

#include "macros_2.h"

double ellK(double, double);


int aav(double *wr, double *wth, double dt)
{
  double lamba_o, F_th, y_th;   // Initial Mino time (for w_th),Incomplete elliptic integral of the firs kind and its argument
  double beta, J;               // Auxiliar variables
  
// Computing W_th(chi)
  beta  = SPIN2*(1.-SQR(ENERGY));
  
  lamba_o = Kcom_th - F_th;
  lamba_o/= sqrt(beta*Z_plus);
  
  if(CHI_NEW >= 0)&&(CHI_NEW <= 0.5*PI)
  {
    y_th = 0.5*PI - CHI_NEW;
    
    F_th = ellK(y_th,KTH);
    lamba_o = Kcom_th - F_th;
    lamba_o/= sqrt(beta*Z_plus);
  
    *w_th = OMEGA_th*lamba_o;
  }
  else if(CHI_NEW >= 0.5*PI)&&(CHI_NEW <= PI)
  {
    y_th = CHI_NEW - 0.5*PI;
    
    F_th = ellK(y_th,KTH);
    lamba_o = Kcom_th - F_th;
    lamba_o/= sqrt(beta*Z_plus);
    
    *w_th = PI - OMEGA_th*lamba_o;
  }
  else if(CHI_NEW >= PI)&&(CHI_NEW <= 1.5*PI)
  {
    y_th = 1.5*PI - CHI_NEW;
    
    F_th = ellK(y_th,KTH);
    lamba_o = Kcom_th - F_th;
    lamba_o/= sqrt(beta*Z_plus);
    
    *w_th = PI + OMEGA_th*lamba_o;
  }
  else if(CHI_NEW >= 1.5*PI)&&(CHI_NEW <= 2.*PI)
  {
    y_th = CHI_NEW - 1.5*PI;
    
    F_th = ellK(y_th,KTH);
    lamba_o = Kcom_th - F_th;
    lamba_o/= sqrt(beta*Z_plus);
    
    *w_th = 2.*PI - OMEGA_th*lamba_o;
  }
  else
  {
    printf("wth geo out of bound");
    exit(0);
  }
  
// Computing W_r(psi)
  
  double wr_psi_save, wr_psi_scale;
  double dwrdpsi, dpsi;
  double r, psi,costh;
  
  wr_psi_save = *wr;
  
  psi = PSI_NEW;
  
  dpsi = PSI_NEW -PSI_OLD;
  
  for (n=1;n<=MAX_NUM_STEPS;n++)
  {
    
    J = (1.-SQR(ENERGY))*(1.-SQR(ECCENTRICITY)) + 2.*(1. - SQR(ENERGY) - (1. -SQR(ECCENTRICITY)/P_p) )*(1.+ECCENTRICITY*cos(psi))
       +( (1.-SQR(ENERGY))*(3.+SQR(ECCENTRICITY))/(1.-SQR(ECCENTRICITY)) -4./P_p
       +(1.-SQR(ECCENTRICITY))/SQR(P_p)*( SPIN2*(1.-SQR(ENERGY)) + SQR(LZ) + Q_CONSTANT))*SQR(1. +ECCENTRICITY*cos(psi));

    dwrdpsi = (1.-SQR(ECCENTRICITY))*OMEGA_r;
    dwrdpsi /= P_p*sqrt(J);
  
    wr_psi_scale = fabs(wr_psi_save) + fabs(dwrdpsi*dpsi) + TINY;
    
    analytic_r_costh(wr_psi_scale, w_th, &r, &costh);
    
    psi = acos((P_p/r - 1.)/ECCENTRICITY);
    
    
    
  }
  
/*----- This is the end -----*/
  return 0;
}

int bs_wr_chi(double *wr_psi_save, double *wTh, double *psi, double *chi, double *phi, double *time, double dt_try, double tolerance,
                       double wR_scale, double wTh_scale, double phi_scale, double *dt_used, double *dt_next)
{
  double wR_save,wTh_save,phi_save;                         /* Variables to save wR_p, wTh_p, and Phi_p               */
  double wR_estimated,wTh_estimated,phi_estimated;          /* Values of wR, wTh, and phi from the midpoint method    */
  double wR_extrapolated,wTh_extrapolated,phi_extrapolated; /* Values of wR, wTh, and phi after extrapolation         */
  double wR_err,wTh_err,phi_err;                            /* Estimated errors from extrapolation                     */
  double tolerance_1;                                        /* Tolerances                                              */
  static double old_tolerance=-1.0;                          /*                                                         */
  double dt;                                                 /* Step                                                    */
  double z_ind;                                              /* Independent variable for polynomial extrapolation       */
  double error_max;                                          /* Normalized error estimate                               */
  double factor,reduction_factor,scale;
  double work,work_min;
  static double time_new;
  
  double error[KSEC];                                        /* Normalized error estimate                               */
  static double a[KSEC];                                     /* Work coefficients                                       */
  static double alpha[KSEC][KSEC];                           /* alpha(k,q)                                              */
  
  static long unsigned k_max,k_optimal;                      /* Optimal row number for convergence                      */
  long unsigned k,kk,iq,km;
  long unsigned reduction;
  static long unsigned first=1;
  long unsigned exit_flag=0;
  
  wR_estimated = 0.0;
  wTh_estimated = 0.0;
  phi_estimated = 0.0;
  
  if (tolerance != old_tolerance)
  {
    *dt_next = time_new = -1.0e29;
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
  
  dt = dt_try;
  
  wR_save = *wR;
  wTh_save = *wTh;
  phi_save = *phi;
  
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
      
      midpoint_method_AAV(k,dt,wR_save,wTh_save,&psi, &chi,phi_save,&wR_estimated,&wTh_estimated,&phi_estimated);
      
      z_ind = dt*dt/k/k/4.0;
      
      poly_extrapol_AAV(k,z_ind,wR_estimated,wTh_estimated,phi_estimated,&wR_extrapolated,&wTh_extrapolated,&phi_extrapolated,&wR_err,&wTh_err,&phi_err);
      
      *wR = wR_extrapolated;
      *wTh = wTh_extrapolated;
      *phi = phi_extrapolated;
      
      if (k != 1)
      {
        error_max = TINY;
        error_max = DMAX(error_max,fabs(wR_err/wR_scale));
        error_max = DMAX(error_max,fabs(wTh_err/wTh_scale));
        error_max = DMAX(error_max,fabs(phi_err/phi_scale));
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


int midpoint_metho(long unsigned n_half_steps, double dt_total, double wR, double wTh, double *psi, double *chi,double phi, double *wR_predicted, double *wTh_predicted, double *phi_predicted)
{
  long unsigned n;
  double z_minus_wR,z_minus_wTh,z_minus_phi;
  double dz_minus_wR,dz_minus_wTh,dz_minus_phi;
  double z_plus_wR,z_plus_wTh/*,z_plus_phi*/;
  double dz_plus_wR,dz_plus_wTh,dz_plus_phi;
  double dz_store_wR,dz_store_wTh,dz_store_phi;
  double derivative_wR,derivative_wTh,derivative_phi;
  double dt,dt_2;
  
  
  /*----- Step Size for this iteration -----*/
  dt = dt_total/(2*n_half_steps);
  
  
  /*----- ZEROth Stage -----*/
  z_minus_wR = wR;
  z_minus_wTh = wTh;
  z_minus_phi = phi;
  
  dz_minus_wR = 0.0;
  dz_minus_wTh = 0.0;
  dz_minus_phi = 0.0;
  
  
  /*----- FIRST Stage -----*/
  osculating_RHS(&derivative_wR,&derivative_wTh,&psi,&chi,&derivative_phi,z_minus_wR,z_minus_wTh);
  
  dz_plus_wR = dt*derivative_wR;
  dz_plus_wTh = dt*derivative_wTh;
  dz_plus_phi = dt*derivative_phi;
  
  z_plus_wR = z_minus_wR + dz_plus_wR;
  z_plus_wTh = z_minus_wTh + dz_plus_wTh;
  /*  z_plus_phi = z_minus_phi + dz_plus_phi;       */
  
  
  /*----- GENERAL Step -----*/
  osculating_RHS(&derivative_wR,&derivative_wTh,&psi,&chi,&derivative_phi,z_plus_wR,z_plus_wTh);
  
  dt_2 = 2.0*dt;
  
  for (n=2;n<=2*n_half_steps;n++)
  {
    dz_store_wR = dz_minus_wR + dt_2*derivative_wR;
    dz_store_wTh = dz_minus_wTh + dt_2*derivative_wTh;
    dz_store_phi = dz_minus_phi + dt_2*derivative_phi;
    
    dz_minus_wR = dz_plus_wR;
    dz_minus_wTh = dz_plus_wTh;
    dz_minus_phi = dz_plus_phi;
    
    dz_plus_wR = dz_store_wR;
    dz_plus_wTh = dz_store_wTh;
    dz_plus_phi = dz_store_phi;
    
    z_plus_wR = z_minus_wR + dz_plus_wR;
    z_plus_wTh = z_minus_wTh + dz_plus_wTh;
    /*    z_plus_phi = z_minus_phi + dz_plus_phi; */
    
    osculating_RHS(&derivative_wR,&derivative_wTh,&psi,&chi,&derivative_phi,z_plus_wR,z_plus_wTh);
  }
  
  dz_store_wR = 0.5*(dz_minus_wR+dz_plus_wR+dt*derivative_wR);
  dz_store_wTh = 0.5*(dz_minus_wTh+dz_plus_wTh+dt*derivative_wTh);
  dz_store_phi = 0.5*(dz_minus_phi+dz_plus_phi+dt*derivative_phi);
  
  *wR_predicted = z_minus_wR + dz_store_wR;
  *wTh_predicted = z_minus_wTh + dz_store_wTh;
  *phi_predicted = z_minus_phi + dz_store_phi;
  
  return 0;
}


int poly_extrapol_AAV(long unsigned call_number,double zi,double wR_mid, double wTh_mid,
                      double phi_mid, double *wR_extrapol, double *wTh_extrapol, double *phi_extrapol, double *dwR, double *dwTh, double *dphi)
{
  double c_wR,c_wTh,c_phi;
  double delta_wR,delta_wTh,delta_phi;
  double q_wR,q_wTh,q_phi;
  double delta,factor_1,factor_2;
  long unsigned j;
  
  zz[call_number] = zi;
  
  *wR_extrapol = wR_mid;
  *wTh_extrapol = wTh_mid;
  *phi_extrapol = phi_mid;
  
  *dwR = wR_mid;
  *dwTh = wTh_mid;
  *dphi = phi_mid;
  
  if (call_number==1)
  {
    tableau[1][1] = wR_mid;
    tableau[2][1] = wTh_mid;
    tableau[3][1] = phi_mid;
  }
  else
  {
    c_wR = wR_mid;
    c_wTh = wTh_mid;
    c_phi = phi_mid;
    
    for (j=1;j<call_number;j++)
    {
      delta = 1.0/(zz[call_number-j]-zi);
      factor_1 = zi*delta;
      factor_2 = zz[call_number-j]*delta;
      
      q_wR = tableau[1][j];
      q_wTh = tableau[2][j];
      q_phi = tableau[3][j];
      
      tableau[1][j] = *dwR;
      tableau[2][j] = *dwTh;
      tableau[3][j] = *dphi;
      
      delta_wR = c_wR - q_wR;
      delta_wTh = c_wTh - q_wTh;
      delta_phi = c_phi - q_phi;
      
      *dwR = factor_1*delta_wR;
      *dwTh = factor_1*delta_wTh;
      *dphi = factor_1*delta_phi;
      
      c_wR = factor_2*delta_wR;
      c_wTh = factor_2*delta_wTh;
      c_phi = factor_2*delta_phi;
      
      *wR_extrapol += *dwR;
      *wTh_extrapol += *dwTh;
      *phi_extrapol += *dphi;
    }
    
    tableau[1][call_number] = *dwR;
    tableau[2][call_number] = *dwTh;
    tableau[3][call_number] = *dphi;
  }
  
  
  /*----- This is the end -----*/
  return 0;
  
}


