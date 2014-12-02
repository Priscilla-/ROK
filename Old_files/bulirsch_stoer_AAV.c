/*==============================================================================

  This routine carries out a Bulirsch-Stoer step with monitoring of local
    truncation error to ensure accuracy and adjust stew_thze.                                        
                                                                                       
  This has been adapted from Numerical Recipes in C.                                   
 -------------------------------------------------------------------------------
  Last Update : 29.12.12  by PCM                                                       
==============================================================================*/

#include<stdio.h>
#include<stdlib.h>

#include "global_quantities_kerr.h"
#include "macros_1.h"

#define K_MAXIMUM 8               //  Maximum row number in the extrapolation              
#define KSEC (K_MAXIMUM+2)        //                                                       
#define SAFETY_1 0.25             //  Safety factors                                       
#define SAFETY_2 0.7              //                                                       
#define TINY 1.0e-30              //  Prevents division by zero                            
#define RED_MIN 0.7               //  Minimum factor for step size reduction               
#define RED_MAX 1.0e-5            //  Maximum factor for step size reduction               
#define SCALE_MAX 0.1             //  1/SCALE_MAX = Maximum factor for step increasing     

int midpoint_method(long unsigned, double, double, double, double, double *, double *, double *);
int poly_extrapol(long unsigned, double, double, double, double, double *, double *, double *, double *, double *, double *);

//int RHs_r_theta_phi(double *, double *, double *, double , double );
int osculating_RHS(double *, double *, double *, double , double);

double tableau[3][KSEC],zz[KSEC];

int bulirsch_stoer_step_AAV(double *w_r, double *w_th, double *phi, double *time, double dt_try, double tolerance, double w_r_scale, double w_th_scale, double phi_scale, double *dt_used, double *dt_next)
{
  double w_r_save,w_th_save,phi_save;                         // Variables to save w_r_f, w_th_f, phi_f               
  double w_r_estimated,w_th_estimated,phi_estimated;          //Values of w_r, w_th, and phi from the midpoint method    
  double w_r_extrapolated,w_th_extrapolated,phi_extrapolated; // Values of w_r, w_th, and phi after extrapolation         
  double w_r_err,w_th_err,phi_err;                            // Estimated errors from extrapolation                     
  double tolerance_1;                                         //Tolerances                                              
  static double old_tolerance=-1.0;                           //                                                         
  double dt;                                                  // Step                                                    
  double z_ind;                                               // Independent variable for polynomial extrapolation       
  double error_max;                                           // Normalized error estimate                               
  double factor,reduction_factor,scale;
  double work,work_min;
  static double time_new;

  double error[KSEC];                                        // Normalized error estimate                               
  static double a[KSEC];                                     // Work coefficients                                       
  static double alpha[KSEC][KSEC];                           // alpha(k,q)                                              

  static long unsigned k_max,k_optimal;                      // Optimal row number for convergence                      
  long unsigned k,kk,iq,km;
  long unsigned reduction;
  static long unsigned first=1;
  long unsigned exit_flag=0;

  w_r_estimated = 0.0;
  w_th_estimated = 0.0;
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
  
  w_r_save  = *w_r;
  w_th_save = *w_th;
  phi_save  = *phi;
  
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

      midpoint_method(k,dt,w_r_save,w_th_save,phi_save,&w_r_estimated,&w_th_estimated,&phi_estimated);
	  
      z_ind = dt*dt/k/k/4.0;
	  
      poly_extrapol(k,z_ind,w_r_estimated,w_th_estimated,phi_estimated,&w_r_extrapolated,&w_th_extrapolated,&phi_extrapolated,&w_r_err,&w_th_err,&phi_err);

      *w_r = w_r_extrapolated;
	  *w_th = w_th_extrapolated;
      *phi = phi_extrapolated;

      if (k != 1)
      {
        error_max = TINY;
     	error_max = DMAX(error_max,fabs(w_r_err/w_r_scale));
		error_max = DMAX(error_max,fabs(w_th_err/w_th_scale));
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


int midpoint_method(long unsigned n_half_steps, double dt_total, double w_r, double w_th, double phi, double *w_r_predicted, double *w_th_predicted, double *phi_predicted)
{
  long unsigned n;
  double z_minus_w_r,z_minus_w_th,z_minus_phi;
  double dz_minus_w_r,dz_minus_w_th,dz_minus_phi;
  double z_plus_w_r,z_plus_w_th;//,z_plus_phi;
  double dz_plus_w_r,dz_plus_w_th,dz_plus_phi;
  double dz_store_w_r,dz_store_w_th,dz_store_phi;
  double derivative_w_r,derivative_w_th,derivative_phi;
  double dt,dt_2;


//----- Step Size for this iteration -----
  dt = dt_total/(2*n_half_steps);      


//----- ZEROth Stage -----
  z_minus_w_r = w_r;      
  z_minus_w_th = w_th;      
  z_minus_phi = phi;      
  
  dz_minus_w_r = 0.0;     
  dz_minus_w_th = 0.0;     
  dz_minus_phi = 0.0;     


//----- FIRST Stage -----
  osculating_RHS(&derivative_w_r,&derivative_w_th,&derivative_phi,z_minus_w_r,z_minus_w_th);
  
  dz_plus_w_r = dt*derivative_w_r;                                                           
  dz_plus_w_th = dt*derivative_w_th;
  dz_plus_phi = dt*derivative_phi;                                                     
        
  z_plus_w_r = z_minus_w_r + dz_plus_w_r;                                                 
  z_plus_w_th = z_minus_w_th + dz_plus_w_th;
//  z_plus_phi = z_minus_phi + dz_plus_phi;       


//----- GENERAL Step -----
  osculating_RHS(&derivative_w_r,&derivative_w_th,&derivative_phi,z_plus_w_r,z_plus_w_th);
  
  dt_2 = 2.0*dt;                                                                              
  
  for (n=2;n<=2*n_half_steps;n++)                                                             
  {
    dz_store_w_r = dz_minus_w_r + dt_2*derivative_w_r;
    dz_store_w_th = dz_minus_w_th + dt_2*derivative_w_th;
    dz_store_phi = dz_minus_phi + dt_2*derivative_phi;

    dz_minus_w_r = dz_plus_w_r;
	dz_minus_w_th = dz_plus_w_th;
    dz_minus_phi = dz_plus_phi;

    dz_plus_w_r = dz_store_w_r;
	dz_plus_w_th = dz_store_w_th;
    dz_plus_phi = dz_store_phi;
	
    z_plus_w_r = z_minus_w_r + dz_plus_w_r;
	z_plus_w_th = z_minus_w_th + dz_plus_w_th;
//    z_plus_phi = z_minus_phi + dz_plus_phi; 

   osculating_RHS(&derivative_w_r,&derivative_w_th,&derivative_phi,z_plus_w_r,z_plus_w_th);
  }

  dz_store_w_r = 0.5*(dz_minus_w_r+dz_plus_w_r+dt*derivative_w_r);
  dz_store_w_th = 0.5*(dz_minus_w_th+dz_plus_w_th+dt*derivative_w_th);
  dz_store_phi = 0.5*(dz_minus_phi+dz_plus_phi+dt*derivative_phi);

  *w_r_predicted = z_minus_w_r + dz_store_w_r;
  *w_th_predicted = z_minus_w_th + dz_store_w_th;
  *phi_predicted = z_minus_phi + dz_store_phi;

  return 0;
}


int poly_extrapol(long unsigned call_number,double zi,double w_r_mid, double w_th_mid, double phi_mid, double *w_r_extrapol, double *w_th_extrapol, double *phi_extrapol, double *dw_r, double *dw_th, double *dphi)
{
  double c_w_r,c_w_th,c_phi;
  double delta_w_r,delta_w_th,delta_phi;
  double q_w_r,q_w_th,q_phi;
  double delta,factor_1,factor_2;
  long unsigned j;

  zz[call_number] = zi;

  *w_r_extrapol = w_r_mid;
  *w_th_extrapol = w_th_mid;
  *phi_extrapol = phi_mid;
  
  *dw_r = w_r_mid;
  *dw_th = w_th_mid;
  *dphi = phi_mid;

  if (call_number==1)
  {
    tableau[1][1] = w_r_mid;
	tableau[2][1] = w_th_mid;
    tableau[3][1] = phi_mid;
  }
  else
  {
    c_w_r = w_r_mid;
	c_w_th = w_th_mid;
    c_phi = phi_mid;
	
    for (j=1;j<call_number;j++)
    {
      delta = 1.0/(zz[call_number-j]-zi);
      factor_1 = zi*delta;
      factor_2 = zz[call_number-j]*delta;

      q_w_r = tableau[1][j];
      q_w_th = tableau[2][j];
	  q_phi = tableau[3][j];
	  
      tableau[1][j] = *dw_r;
      tableau[2][j] = *dw_th;
	  tableau[3][j] = *dphi;
	  
      delta_w_r = c_w_r - q_w_r;
	  delta_w_th = c_w_th - q_w_th;
      delta_phi = c_phi - q_phi;
	  
      *dw_r = factor_1*delta_w_r;
	  *dw_th = factor_1*delta_w_th;
      *dphi = factor_1*delta_phi;
	  
      c_w_r = factor_2*delta_w_r;
	  c_w_th = factor_2*delta_w_th;
      c_phi = factor_2*delta_phi;
	  
      *w_r_extrapol += *dw_r;
	  *w_th_extrapol += *dw_th;
      *phi_extrapol += *dphi;
    }
	
    tableau[1][call_number] = *dw_r;
	tableau[2][call_number] = *dw_th;
    tableau[2][call_number] = *dphi;
  }
  
  
//----- This is the end -----  
  return 0;
  
}
