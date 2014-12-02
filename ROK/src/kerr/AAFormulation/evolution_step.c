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
#include <math.h> 

#include "global_quantities_kerr.h"
#include "global_quantities_main.h"
#include "global_quantities_osc.h"
#include "physical_quantities.h"
#include "macros_2.h"
//#include "macros_1.h"

#define TINY 1.0e-30
#define MAX_NUM_STEPS 10000000


double *Realvector(const long nl, const long nh);
void free_Realvector(double *v, const long nl, const long nh);

int AA_RHS(double *,double*);
int geodesic_evolution_AV(double *x, double *dxdt);
int bulirsch_stoer_AAV(double *, double *, double, double, double *, double *, double *);

int evolution_step(double *x, double t_initial, double t_final,  double dt_min, double dt_initial,double accuracy)
{
  double *var_save;
  double *var_scale;       // Scaling of the errors
  double dt_used,dt_next, dt, time;  // Time step used by the BS routine and the one computed for the next step
  double *dxdt;            // vector of evolving variables
  long unsigned n;
  int i;
  
  var_scale = Realvector(1, NVAR);
  var_save  = Realvector(1, NVAR);
  dxdt = Realvector(1, NVAR);
    
  time = t_initial;
  dt = SIGN(dt_initial,t_final-t_initial);

  for (i = 1; i <= NVAR; i++)
    var_save[i] = x[i];
 
    
  for (n=1;n<=MAX_NUM_STEPS;n++)
  {
     if(MAG_SF == 0)
        geodesic_evolution_AV(var_save, dxdt);
    else
        AA_RHS(var_save,dxdt);
   
    for (i = 1; i <= NVAR; i++)
      var_scale[i] = fabs(var_save[i]) + fabs(dxdt[i]*dt) + TINY;
     
      if ( (time+dt-t_final)*(time+dt-t_initial) > 0.0 ) // Decreases the time step if it is going to overshoot 
       dt = t_final - time;

//--Calling the integrator
    bulirsch_stoer_AAV(var_save,&time,dt,accuracy,var_scale,&dt_used,&dt_next);
    
    // printf("var_scale[6]=%4.6e var_scale[7] =%4.6e\n",var_scale[6],var_scale[7]);
    //printf("var_save[6]=%4.6e var_save[7] =%4.6e\n",var_save[6],var_save[7]);
     
    if(ODEcount == 1)
    return 0;
    
    if ( (time-t_final)*(t_final-t_initial) >= 0)
    {
      for (i = 1; i <= NVAR; i++)
        x[i] = var_save[i];
        
       free_Realvector(var_save, 1, NVAR);
       free_Realvector(var_scale, 1, NVAR);
       free_Realvector(dxdt,1,NVAR);
        
      return 0;
    }
    if (fabs(dt_next) <= dt_min)
    {
      printf("  | Step size too small for computing derivatives!!\n");
      ODEcount = 1;
      return 0;
    }    
    dt = dt_next;
  }
    
  ODEcount = 1;
  printf("  | Too many steps!!!\n");
    
/*----- This is the end -----*/
    return 0;
}






