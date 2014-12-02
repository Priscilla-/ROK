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

#include "global_quantities_kerr.h"
#include "global_quantities_main.h"

#include "physical_quantities.h"
#include "macros_2.h"

#include "macros_1.h"

#include "global_quantities_osc.h"
#include "physical_quantities.h"
#include "macros_2.h"

#define K_MAXIMUM 9        // Maximum row number in the extrapolation
#define KSEC (K_MAXIMUM+2)  //
#define SAFETY_1 0.25       // Safety factors
#define SAFETY_2 0.7        //
#define TINY 1.0e-30        // Prevents division by zero
#define RED_MIN 0.7         // Minimum factor for step size reduction
#define RED_MAX 1.0e-5      // Maximum factor for step size reduction
#define SCALE_MAX 0.1       // 1/SCALE_MAX = Maximum factor for step increasing


double tableau[7][KSEC],zz[KSEC];

int geodesic_evolution_AV(double *, double*);
int AA_RHS(double *,double*);

int midpoint_method_AAV(long unsigned n_half_steps, double dt_total, double *x, double *x_predicted);
int poly_extrapol_AAV(long unsigned call_number,double zi,double *x_mid, double *x_extrapol, double *dx);


double *Realvector(const long nl, const long nh);
void free_Realvector(double *v, const long nl, const long nh);

int bulirsch_stoer_AAV(double *x, double *time, double dt_try, double tolerance, double *x_scale, double *dt_used, double *dt_next)
{
   int i;
   double *x_save;             // Variables to save Chi_p, Psi_p, and Phi_p
   double *x_estimated;   // Values of chi, psi, and phi from the midpoint method
   double *x_extrapolated; // Values of chi, psi, and phi after extrapolation
   double *x_err; // Estimated errors from extrapolation
  
   x_save = Realvector(1, NVAR);
   x_estimated = Realvector(1, NVAR);
   x_extrapolated = Realvector(1, NVAR);
   x_err = Realvector(1, NVAR);
  
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
  
   for (i = 1; i <= NVAR; i++)
    x_estimated[i] = 0.0;
    
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
      alpha[k][iq] = pow(tolerance_1,(a[k+1] - a[iq+1])/((a[iq+1] - a[1] + 1.0)*(2*k+1)));
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
 
    
   for (i = 1; i <= NVAR; i++)
    x_save[i] = x[i];
  
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
      
      midpoint_method_AAV(k,dt,x_save,x_estimated);
      
      z_ind = dt*dt/k/k/4.0;
      
      poly_extrapol_AAV(k,z_ind,x_estimated,x_extrapolated,x_err);
      
      for (i = 1; i <= NVAR; i++)
       x[i] = x_extrapolated[i];
      
      if (k != 1)
      {
      error_max = TINY;
      for (i = 1; i <= NVAR; i++)
       error_max = DMAX(error_max,fabs(x_err[i]/x_scale[i]));
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
  
  
  free_Realvector(x_err, 1, NVAR);
  free_Realvector(x_estimated, 1, NVAR);
  free_Realvector(x_save, 1, NVAR);
  free_Realvector(x_extrapolated, 1, NVAR);

   return 0;
}


int midpoint_method_AAV(long unsigned n_half_steps, double dt_total, double *x, double *x_predicted)
{
  long unsigned n;
  double *z_minus_x, *dz_minus_x, *z_plus_x, *dz_plus_x,*dz_store_x;
  double dt,dt_2;
  int i;
  
  z_minus_x = Realvector(1, NVAR);
  dz_minus_x = Realvector(1, NVAR);
  z_plus_x = Realvector(1, NVAR);
  dz_plus_x = Realvector(1, NVAR);
  dz_store_x = Realvector(1, NVAR);
   
/*----- Step Size for this iteration -----*/
   dt = dt_total/(2*n_half_steps);
   
/*----- ZEROth Stage -----*/
  for (i = 1; i <= NVAR; i++)
  {
    z_minus_x[i] = x[i];
    dz_minus_x[i] = 0.0;
  }
  
  //  printf("z_minus_x[6]=%4.6e z_minus_x[7] =%4.6e\n",z_minus_x[6],z_minus_x[7]);
  //  printf("x_predicted[6]=%4.6e x_predicted[7] =%4.6e\n",x_predicted[6],x_predicted[7]);

/*----- FIRST Stage -----*/
  if(MAG_SF == 0)
    geodesic_evolution_AV(z_minus_x,x_predicted);
  else
    AA_RHS(z_minus_x,x_predicted);
 
  for (i = 1; i <= NVAR; i++)
  {
    dz_plus_x[i] = dt*x_predicted[i];
    z_plus_x[i]  = z_minus_x[i] + dz_plus_x[i];
    
   }
  // printf("z_plus_x[6]=%4.6e z_plus_x[7] =%4.6e\n",z_plus_x[6],z_plus_x[7]);
  // printf("x_predicted[6]=%4.6e x_predicted[7] =%4.6e\n",x_predicted[6],x_predicted[7]);

/*----- GENERAL Step -----*/
  if(MAG_SF == 0)
    geodesic_evolution_AV(z_plus_x,x_predicted);
  else
    AA_RHS(z_plus_x,x_predicted);

  dt_2 = 2.0*dt;
  
  // printf("dt_2 =%4.6e\n",dt_2);
  
  for (n=2;n<=2*n_half_steps;n++)
  {
    for (i = 1; i <= NVAR; i++)
    {
     dz_store_x[i] = dz_minus_x[i] + dt_2*x_predicted[i];
     dz_minus_x[i] = dz_plus_x[i];
     dz_plus_x[i]  = dz_store_x[i];
     z_plus_x[i]   = z_minus_x[i] + dz_plus_x[i];
    }
    
  /*   printf(" x_predicted[6]=%4.6e  x_predicted[7] =%4.6e\n", x_predicted[6], x_predicted[7]);
     printf(" dz_store_x[6]=%4.6e  dz_store_x[7] =%4.6e\n", dz_store_x[6], dz_store_x[7]);
     printf(" dz_minus_x[6]=%4.6e  dz_minus_x[7] =%4.6e\n", dz_minus_x[6], dz_minus_x[7]);
     printf(" dz_plus_x[6]=%4.6e   dz_plus_x[7] =%4.6e\n", dz_plus_x[6], dz_plus_x[7]);
    
    printf("z_plus_x[6]=%4.6e z_plus_x[7] =%4.6e\n",z_plus_x[6],z_plus_x[7]);
    printf("x_predicted[6]=%4.6e x_predicted[7] =%4.6e\n",x_predicted[6],x_predicted[7]);
*/
    if(MAG_SF == 0)
      geodesic_evolution_AV(z_plus_x,x_predicted);
    else
      AA_RHS(z_plus_x,x_predicted);
    
  }
    
  for (i = 1; i <= NVAR; i++)
  {
    dz_store_x[i] = 0.5*(dz_minus_x[i] + dz_plus_x[i] + dt*x_predicted[i]);
    x_predicted[i] = z_minus_x[i] + dz_store_x[i];
  }
  
  free_Realvector(z_minus_x, 1, NVAR);
  free_Realvector(dz_minus_x, 1, NVAR);
  
  free_Realvector(z_plus_x, 1, NVAR);
  free_Realvector(dz_plus_x, 1, NVAR);
  
  free_Realvector(dz_store_x, 1, NVAR);
 
  return 0;
}


int poly_extrapol_AAV(long unsigned call_number,double zi,double *x_mid, double *x_extrapol, double *dx)
{
   double *c_x, *delta_x, *q_x;
   double delta,factor_1,factor_2;
   long unsigned j;
   int i;
  
   c_x     = Realvector(1, NVAR);
   delta_x = Realvector(1, NVAR);
   q_x     = Realvector(1, NVAR);
  
   zz[call_number] = zi;
  
  for (i = 1; i <= NVAR; i++)
  {
    x_extrapol[i] = x_mid[i];
    dx[i] = x_mid[i];
   
   if (call_number==1)
   {
     tableau[i][1] = x_mid[i];
     
   }
   else
   {
     c_x[i] = x_mid[i];
     
     for (j=1;j<call_number;j++)
     {
      delta = 1.0/(zz[call_number-j]-zi);
      factor_1 = zi*delta;
      factor_2 = zz[call_number-j]*delta;
      
      q_x[i] = tableau[i][j];
      tableau[i][j] = dx[i];
      delta_x[i] = c_x[i] - q_x[i];
       
      dx[i] = factor_1*delta_x[i];
       
      c_x[i] = factor_2*delta_x[i];
       
      x_extrapol[i] += dx[i];
     }
     
     tableau[i][call_number] =  dx[i];
   }
  }
  
  free_Realvector(c_x, 1, NVAR);
  free_Realvector(q_x, 1, NVAR);
  free_Realvector(delta_x, 1, NVAR);

   return 0;
}
