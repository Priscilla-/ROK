/*==============================================================================
 
 This routine carries out a Bulirsch-Stoer step with monitoring of local
 truncation error to ensure accuracy and adjust stewThze.
 
 This has been adapted from Numerical Recipes in C.
 -------------------------------------------------------------------------------
 Last Update : 17.08.13  by PCM
 ==============================================================================*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>


#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"
#include "macros_1.h"
#include "macros_2.h" 

#define K_MAXIMUM 8              /*  Maximum row number in the extrapolation              */
#define KSEC (K_MAXIMUM+2)        /*                                                       */
#define IMAXX (K_MAXIMUM + 1)
#define SAFETY_1 0.25             /*  Safety factors                                       */
#define SAFETY_2 0.7              /*                                                       */
#define TINY 1.0e-30              /*  Prevents division by zero                            */
#define RED_MIN 0.7               /*  Minimum factor for step size reduction               */
#define RED_MAX 1.0e-5            /*  Maximum factor for step size reduction               */
#define SCALE_MAX 0.1            // 1/SCALE_MAX = Maximum factor for step increasing
#define TINY 1.0e-30             // Prevents division by zero

int checks_main_errors(int);
int osculating_RHS(double *,double*);

double *Realvector(const long nl, const long nh);
double **Realmatrix(const long nrl, const long nrh, const long ncl, const long nch);
void free_Realvector(double *v, const long nl, const long nh);
void free_Realmatrix(double **m,const long nrl, const long nrh,const long ncl, const long nch);




int midpoint_method_AAV(double *x, double *dxdt,double *time, double dt, int nstep, double *z_derivative)
{
    int n,i;
    double *z_minus,*z_plus,z_store;
    double t_total,dt_2,h;
    
    z_minus = Realvector(1, NVAR);
    z_plus  = Realvector(1, NVAR);
    
    h = dt/nstep;
    printf("time= %4.6e dt=%4.6e nstep =%d\n", *time, dt,nstep);

    for (i = 1; i <= NVAR; i++)
    {
        z_minus[i] = x[i];
        z_plus[i]  = x[i] + h*dxdt[i];
    }
    t_total = *time + h;
    
    
    osculating_RHS(z_derivative,z_plus);
    
    dt_2 = 2.0*h;
    
    for (n=2;n<=nstep;n++)
    {
        for (i = 1; i <= NVAR; i++)
        {
            z_store    = z_minus[i] + dt_2*z_derivative[i];
            z_minus[i] = z_plus[i];
            z_plus[i]  = z_store;
        }
        t_total += h;
        
        osculating_RHS(z_derivative,z_plus);
    }
    for (i = 1; i<= NVAR; i++)
        z_derivative[i] = 0.5*(z_minus[i] + z_plus[i] + h*z_derivative[i]);
    
    free_Realvector(z_minus, 1, NVAR);
    free_Realvector(z_plus, 1, NVAR);
    
    return 0;
}


int poly_extrapol_AAV(long unsigned call_number,double zi,double *z_mid,double *z_extrapol, double *dz, double **tableau, double *err)
{
    double *fx;
    double xx,q_z,q_z1,c;
    double delta,ddz;
    long unsigned j;
    int i;
    
    fx = Realvector(1, call_number);
    
    err[call_number] = zi;
    
    if (call_number==1)
    {
        for(i=1; i<=NVAR; i++)
        {
            z_extrapol[i] = z_mid[i];
            tableau[i][1] = z_mid[i];
            dz[i] = z_mid[i];
        }
    }
    else
    {
        for (j=1;j<call_number;j++)
            fx[j+1] = err[call_number-j]/zi;
        
        for(i=1; i<=NVAR; i++)
        {
            delta = tableau[i][1];
            
            c = z_mid[i];
            xx = c;
            tableau[i][1] = xx;
            
            for (j=2;j<=call_number;j++)
            {
                q_z1 = fx[j]*delta;
                q_z  = q_z1 - c;
                
                if(q_z !=0.)
                {
                    q_z = (c-delta)/q_z;
                    ddz = c*q_z;
                    c = q_z1*q_z;
                }
                else
                    ddz = delta;
                
                if (j != call_number)
                    delta = tableau[i][j];
                
                tableau[i][j] = ddz;
                xx += ddz;
            }
            dz[i] = ddz;
            z_extrapol[i] = xx;
        }
    }
    
    free_Realvector(fx,1, call_number);
    
    return 0;
}



int bulirsch_stoer_AAV(double *x,double *dxdt, double *time, double dt_try, double tolerance, double *x_scale, double *dt_used, double *dt_next)
{
  
  double *var_estimated; // Values used in the midpoint method
  double *var_error;     // Errors estimated from extrapolation
  double *var_save;
  double *var_x;
  double tolerance_1;                                        /* Tolerances                                              */
  double dt;                                                 /* Step                                                    */
  double z_ind;                                              /* Independent variable for polynomial extrapolation       */
  double error_max;                                          /* Normalized error estimate                               */
  double factor,reduction_factor,scale;
  double work,work_min;
  double **tableau,*zz;
  
  static double old_tolerance=-1.0;                          /*                                                         */
  static double time_new;
  static double a[IMAXX +1];                                     /* Work coefficients                                       */
  static double alpha[K_MAXIMUM +1][K_MAXIMUM +1];                           /* alpha(k,q)                                              */
  static int nseq[IMAXX +1] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18};
  static int k_max,k_optimal;                      /* Optimal row number for convergence                      */
  static int first=1;


  int  exit_flag=0;
  int  k,kk,iq,km;
  int  reduction,i;
  
  tableau = Realmatrix(1, NVAR, 1, K_MAXIMUM);
  zz = Realvector(1, K_MAXIMUM);
  
  var_x = Realvector(1,K_MAXIMUM);
  var_error = Realvector(1, NVAR);
  var_save = Realvector(1, NVAR);
  var_estimated = Realvector(1, NVAR);
  
  if (tolerance != old_tolerance)
  {
    *dt_next = time_new = -1.0e29;
    tolerance_1 = SAFETY_1*tolerance;
    a[1] = nseq[1] + 1;
    
    for (k=1;k<=K_MAXIMUM;k++)
    { a[k + 1] = a[k] + nseq[k + 1]; printf("a = %4.6e\n", a[k+1]);}
      
      
    
    for (iq=2;iq<=K_MAXIMUM;iq++)
    {
      for (k=1;k<iq;k++)
      { alpha[k][iq] = pow(tolerance_1,(a[k+1] - a[iq+1])/((a[iq+1] - a[1] + 1.0)*(2*k+1)));
          
         printf("alpha = %4.6e tol = %4.6e\n", alpha[k][iq],tolerance_1);      }
    }
    old_tolerance = tolerance;
    
    for (k_optimal=2;k_optimal<K_MAXIMUM;k_optimal++)
     if (a[k_optimal+1]>a[k_optimal]*alpha[k_optimal-1][k_optimal]) break;
    k_max = k_optimal;
  }
  dt = dt_try;

printf("dt=%4.6e old_tolerance=%4.6e\n",dt,old_tolerance);

  for (i = 1; i <= NVAR; i++)
   var_save[i] = x[i];
  
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
      printf("time_new= %4.6e k =%d\n", time_new, k);
        if (time_new == (*time))
        {
            ODEcount =1;
            printf("  | Step size underflow!!!\n");exit(0);
        }
        printf("time= %4.6e dt=%4.6e nstep =%d\n", *time, dt,nseq[k]);

      midpoint_method_AAV(var_save,dxdt,*time, dt, nseq[k], var_estimated);
        
      
      z_ind =SQR(dt/nseq[k]);
        
        
        printf("nseq[k]= %d dt=%4.6e z_ind =%4.6e\n", nseq[k], dt,z_ind);
        
        
        
      poly_extrapol_AAV(k,z_ind,var_estimated,x,var_error, tableau, var_x);
      
      if (k != 1)
      {
        error_max = TINY;
        for (i=1;i<=NVAR;i++)
         error_max = DMAX(error_max,fabs(var_error[i]/x_scale[i]));
         error_max /= tolerance;
         km  = k - 1;
         zz[km] = pow(error_max/SAFETY_1,1.0/(2*km+1));
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
            reduction_factor = SAFETY_2/zz[km];
            break;
          }
          else if (k == k_optimal && alpha[k_optimal-1][k_optimal]<zz[km])
          {
             reduction_factor = 1.0/zz[km];
             break;
          }
          else if (k_optimal == k_max && alpha[km][k_max-1]<zz[km])
          {
             reduction_factor = alpha[km][k_max-1]*SAFETY_2/zz[km];
              break;
          }
          else if (alpha[km][k_optimal]<zz[km])
          {
              reduction_factor = alpha[km][k_optimal-1]/zz[km];
              break;
           }
         }
        }
      
      printf("redfac =%4.6e\n", reduction_factor);
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
        factor = DMAX(zz[kk],SCALE_MAX);
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
  
  free_Realvector(var_error, 1, NVAR);
  free_Realvector(var_estimated, 1, NVAR);
  free_Realvector(var_save, 1, NVAR);
  free_Realvector(var_x, 1, K_MAXIMUM);
  
  free_Realvector(zz, 1, K_MAXIMUM);
  free_Realmatrix(tableau, 1, K_MAXIMUM, 1, K_MAXIMUM);

  return 0;
}


