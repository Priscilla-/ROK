/*==============================================================================
 
 This routine evolves the geodesic equations for the Kerr Spacetime using 
 osculating evolution conditions (see my notes) together with an action angle 
 variable formulation.
 
 NOTE: Here we are effectively computing the:

 Dissipative SF: Gives the secular change on E, Lz and Q. Mimics the RR effects 
 by its time evolution.
 
 Conservative SF: Shiftes the instantaneous values of the orbital parameters and
 changes the orbital phases
 
 -------------------------------------------------------------------------------
                                            Created by Priscilla on 20/12/2012
 -------------------------------------------------------------------------------
 Last Update: 31.07.13
 =============================================================================*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/types.h>
#include<time.h>

#include "global_quantities_kerr.h"
#include "global_quantities_main.h"
#include "global_quantities_osc.h"
#include "physical_quantities.h"
#include "macros_2.h" 

int analytic_r_costh(double*,double *,double *);
int check_if_orbit_plunges(double, double);
int compute_frequencies(double *, double *, double *, double *, double *,double *, double *, double *, double *);
int compute_Keplerian_period(double, double, double *, long unsigned *, double);
int evolution_step(double *,double, double, double, double, double );
int geodesic_evolution(double *, double *, double *,double *, double, double, double, double, double);
//int initialization(void);

double *Realvector(const long nl, const long nh);
void free_Realvector(double *v, const long nl, const long nh);

int evolution_OS(FILE *data, FILE *dataf, FILE *dataflux, double *x)
{
   double step_insp,StepsKTsec;
 // double omegaR_OLD, omegaTH_OLD, omegaT_OLD;
   double dt_sec, t_evol,M;                                    // vector of variables
   int i;
  long unsigned N_wr, N_wth; // Number of turns
  
//----Initial fundamental frequencies
//  omegaR_OLD  = OMEGA_r;
//  omegaTH_OLD = OMEGA_th;
//  omegaT_OLD  = OMEGA_t;
  
//----Time interval for data saving
//  dt_sec = ceil(STORE_rate*Nfreq );
  
    double dimensional_factor = MASS_SUN*GN/(LIGHT_SPEED*SQR(LIGHT_SPEED));
    
    dt_sec = DTsec;//STORE_rate;
  step_insp = (int) (TIME_yrs*SECSYR/dt_sec);
  if (step_insp < 3) step_insp = 3;

/*----- Keplerian period and number of time steps inside it-----*/
  if(MAG_SF == 0)
    compute_Keplerian_period(P_p,ECCENTRICITY,&KEPLER_period,&StepKT, dt_sec); 

/*----MBH mass in sec----*/
  M = MASS_in_MSUN*MASS_SUN*(GN/CSPEED3); 
  KEPLER_period *= sqrt(M)*M;
  StepsKTsec  = KEPLER_period/dt_sec  ;

 // printf("Insp = %4.6e [sec], dt = %4.6e [sec], #stepInsp = %4.6e, Period (Keplerian) = %4.6e [sec], #stepKep = %4.6e\n ", TIME_yrs*SECSYR, dt_sec, step_insp, KEPLER_period,StepsKTsec);
  
  t_evol = 0.;
   
      
  for (i = 0; i <= step_insp; i++)
  {

      
    if(MAG_SF != 0)
    {
      PSI_NEW = PSI_OLD;
      CHI_NEW = CHI_OLD;
      x[1] = PHI_OLD;
      x[2] = TAU;
      x[3] = P_p;
      x[4] = ECCENTRICITY;
      x[5] = INCLINATION;
      x[6] = Qr;
      x[7] = Qth;
    }
    else
    {
      x[1] = Qr;
      x[2] = Qth;
      x[3] = PHI_OLD;
      x[4] = TIME;
    }
    
    
    analytic_r_costh(x, &R_p_o, &COSTH_p_o);
    
      
    if (STORE_rate!= 0)
    {
      if(MAG_SF != 0)
      {
        fprintf(dataflux,"%4.6e  %4.12e  %4.12e  %4.14e\n",t_evol, ENERGY_dot,LZ_dot,Q_CONSTANT_dot);
        fprintf(dataf,"%4.6e %4.6e %4.6e %4.6e %4.6e %4.6e %4.6e %4.6e %4.6e\n",t_evol, x[1], x[2], x[3], x[4], x[5], x[6], x[7]);
        fprintf(data, "%4.6e %4.6e %4.6e %4.6e \n",t_evol, R_p_o, COSTH_p_o, x[3]);
      }
      else
      {
        fprintf(dataf,"%4.6e %4.6e %4.6e %4.6e \n",t_evol, x[1], x[2], x[3]);
        fprintf(data, "%4.6e %4.6e %4.6e %4.6e \n",t_evol, R_p_o, COSTH_p_o, x[3]);
      }
    }
    
     printf("%4.6e %4.6e %4.6e %4.6e \n",t_evol, R_p_o, COSTH_p_o , x[3]);
    evolution_step(x, t_evol, t_evol+dt_sec,DT_MIN/BHMASS_to_s,dt_sec, ACCURACY/BHMASS_to_s);
     
    /**----- Check whether the system has plunged -----**/
    check_if_orbit_plunges(ECCENTRICITY,Q_CONSTANT);
    

    if(MAG_SF != 0)
    {
     PSI_OLD = PSI_NEW;
     CHI_OLD = CHI_NEW;
     PHI_OLD = x[1];
     TAU = x[2];
     P_p = x[3];
     ECCENTRICITY = x[4];
     INCLINATION = x[5];
     Qr = x[6];
     Qth = x[7];
    }
    else
    {
     Qr = x[1];
     Qth= x[2];
     PHI_OLD= x[3];
     TIME =x[4];
     }
  
   t_evol += dt_sec;
  }
/*======= EVOLUTION  LOOP  ENDS  HERE! ========*/
    
/*----- This is the end -----*/  
  return 0;
}
