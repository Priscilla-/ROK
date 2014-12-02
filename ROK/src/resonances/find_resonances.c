/*=====================================================================================
 
 For a given resonant order (k,n), such that k*Omega_th + n*Omega_r = 0, and fixed BH
 spin (a/M.), initial eccentricity (e) and inclination (theta_inc), this routine finds the
 (e,p,theta_inc) on resonance,that is for a  given resonant geodesic.
 --------------------------------------------------------------------------------------
                                                    Created by Priscilla on 28/12/2012
 --------------------------------------------------------------------------------------
 Last Update : 06.10.13 by PCM
 ======================================================================================*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"
#include "macros_2.h"

int initialization(void);
int search_resonant_frequencies(double, double, double,  FILE *);

int find_resonances(void)
{
  FILE *data;
  
  long unsigned i,j,k;                  // Counters
  double Ne, Np, Ns;                  // Maximum Numer of iterations in the searching algorithm
  
  double Ec_o_search, P_o_search;     // Initial values of the eccentricity and semilatus rectum in resonance search
  double EC_f_search, P_f_search;     // Final values of the eccentricity and semilatus rectum in resonance search
  double DEC, DP;                     // Eccentricity and semi-latus rectum step size
  double e_search, p_search, a_search;// (e,p,Th_inc) used in a resonance search
  double he, hp;                      // Searching "distance" for e and p
  double theta_inc_search;
  double d_th;
  int Nth;
  
/*----- Initializing some quantities -----*/
  initialization();

  Ec_o_search = ECCENTRICITY;
  P_o_search  = P_p;
  EC_f_search = EMRI_Parameter[15];
  P_f_search  = EMRI_Parameter[16];              
  DEC = Real_Parameter[6];
  DP = Real_Parameter[7];
  
  RES_COUNT = 0;
 
  e_search = Ec_o_search;
       
  he = fabs( EC_f_search - Ec_o_search );
  hp = fabs( P_f_search - P_o_search );
    
  Ne = he/DEC;
  Np = hp/DP;
  
/*--- Search in Theta_inc [positive angles only]---*/
  d_th = 5.0*PI/180;
  Nth  = 18.0;
  theta_inc_search = d_th;
      
/*===== Searching resonances ======*/
   
 data = fopen(FN_Res_frequencies,"a");

 for (k =0; k < Nth; k++)
 {
  if (theta_inc_search < PI/2)
  {
   for(i=0; i<=Ne; i++)
   {
    if(e_search < EC_f_search)
      e_search = EC_f_search;
      
      p_search = P_o_search;
        
     for(j= 0; j<= Np ; j++)
     {
        if(p_search < P_f_search)
         p_search = P_f_search;
         
          search_resonant_frequencies(e_search, p_search, theta_inc_search , data);
  
        p_search -= DP;
      }
      e_search -= DEC;
    }
  }
   theta_inc_search +=d_th;
 }
  
  if( RES_COUNT == 0 )
  {
    printf("\n");
    printf("---------------------------------------------------------------------------\n\n");
    printf("[!] NO RESONANCE has been found within the selected range of parameters :( \n\n");
      
    fprintf(data,"------------------------------------------------------------------------\n\n");
    fprintf(data,"[!] No resonance has been found within the following range of parameters:\n\n");
    fprintf(data,"e_o = %4.4e, e_f = %4.4e, p_o = %4.4e, p_f = %4.4e, th_inc = %4.4e\n\n",
                Ec_o_search, e_search, P_o_search, p_search, theta_inc_search);
    fprintf(data,"------------------------------------------------------------------------\n\n");
  }
    
 fclose(data);

    
 if(EVAL_arg == 1)
 {
    printf("\n");
    printf("In this search we are crossed invalid arguments for the Elliptical integrals\n\n");
   exit(0);
 }

/*----- This is the end -----*/
    return 0;
}
