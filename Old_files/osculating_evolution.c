/*==============================================================================
 This routine Integrates the osculating conditions for the Orbital frequencies 
 in terms of observer's time. 
 -------------------------------------------------------------------------------
                                            Created by Priscilla on 26/03/2013
 -------------------------------------------------------------------------------
 Last Update: 08/04/2013                    Modified by PCM
 =============================================================================*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"

#include "physical_quantities.h"
#include "macros_2.h"

#define TINY 1.0e-30
#define MAX_NUM_STEPS 200000
                    //10000000

int osculating_RHS(double *, double *, double *, double, double);
int bulirsch_stoer_step_AAV(double *, double *, double *, double *, double, double, double, double, double, double *, double *);

int osculating_evolution(double *w_r_o, double *w_th_o, double *phi_o, double t_initial, double t_final, double dt_initial, double dt_min, double accuracy)
{
       
  double w_r_save,w_th_save,phi_save;              // AA Variables to save 
  double w_r_dot, w_th_dot, phi_dot;               // Derivatives of wr, wth, and phi
  double w_r_f, w_th_f, phi_f;                     // Scaling of the errors
  double time,dt ;                                 // Time and Time-step size
  double dt_used,dt_next;                          // Time step used by the BS routine and the one computed for the next step
  long unsigned n;
    
  time = t_initial;
  dt = SIGN(dt_initial,t_final-t_initial);
    
  w_r_save  = *w_r_o;
  w_th_save = *w_th_o;
  phi_save  = *phi_o;
  
  for (n=1;n<=MAX_NUM_STEPS;n++)
  {
    osculating_RHS(&w_r_dot, &w_th_dot, &phi_dot, w_r_save, w_th_save);
     
    w_r_f  = fabs(w_r_save)  + fabs(w_r_dot*dt)  + TINY;
    w_th_f = fabs(w_th_save) + fabs(w_th_dot*dt) + TINY;
    phi_f  = fabs(phi_save)  + fabs(phi_dot*dt)  + TINY;
    
    /*  printf("a)  %4.6e   %4.6e   %4.6e  %4.6e   %4.6e   %4.6e %4.6e % 4.6e   %4.6e %4.6e\n",
             dt , w_r_f, fabs(w_r_save),fabs(w_r_dot*dt),
                 w_th_f, fabs(w_th_save) ,fabs(w_th_dot*dt),
                 phi_f , fabs(phi_save)  , fabs(phi_dot*dt) );
     */
      
      
      if ( (time+dt-t_final)*(time+dt-t_initial) > 0.0 )  /* Decreases the time step if it is going to overshoot */
      {
          dt = t_final - time; printf("HOLA\n");
      }
      
      bulirsch_stoer_step_AAV(&w_r_save,&w_th_save,&phi_save,&time,dt,accuracy,w_r_f,w_th_f,phi_f,&dt_used,&dt_next);
    
      
      if ( (time-t_final)*(t_final-t_initial) >= 0)
      {
        *w_r_o  = w_r_f;
        *w_th_o = w_th_f;
        
        *phi_o  = phi_f; //printf("JOLIN\n");
        
        return 0;
     }
        
      if (fabs(dt_next) <= dt_min){ printf("  | Step size too small for computing derivatives!!\n");
          printf("\n BS %4.6e   %4.6e \n",fabs(dt_next), dt_min);
      }
        
     dt = dt_next;
      
      printf("time=%4.6e dt=%4.6e t_final=%4.6e  t_initial=%4.6e  dt_used=%4.6e  dt_next=%4.6e \n",time, dt, t_final, t_initial,dt_used, dt_next );
    }
    
    printf("  | Too many steps!!!\n");exit(0);
    

/*----- This is the end -----*/
 return 0;
}








