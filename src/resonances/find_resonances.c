/*=====================================================================================
 
 For a given resonant order (k,n), such that k*Omega_th + n*Omega_r = 0, and fixed BH
 spin (a/M.), initial eccentricity (e) and inclination (theta_inc), this routine finds the
 (e,p,theta_inc) on resonance,that is for a  given resonant geodesic.
 --------------------------------------------------------------------------------------
                                                    Created by Priscilla on 28/12/2012
 --------------------------------------------------------------------------------------
 Last Update : 13.07.13 by PCM
 ======================================================================================*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"
#include "macros_2.h"

int initialization(void);
int search_resonant_frequencies(double, double, double, FILE *);

int find_resonances(void)
{
  FILE *data;
  
  long unsigned i,j;              // Counters
  double Ne, Np, Ns;                  // Maximum Numer of iterations in the searching algorithm
  
  double e_search, p_search, a_search;  // (e,p,Th_inc) used in a resonance search
  double he, hp;                                // Searching distance for e and p
    
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
  
 /* d_th = 5.0*PI/180;
  Nth  = 18.0;
  theta_inc_search = d_th;
  theta_inc_search = THETA_inc;

 */
    
  Ns = 18;
  a_search = SPIN;
    
   // printf("a_search=%4.5e\n",a_search);exit(0);
    
/*===== Searching resonances ======*/
   
 data = fopen(FN_Res_frequencies,"a");

  for(i=0; i<=Ne; i++)
  {
    if(e_search < EC_f_search)
      e_search = EC_f_search;
      
      p_search = P_o_search;
        
     for(j= 0; j<= Np ; j++)
     {
        if(p_search < P_f_search)
         p_search = P_f_search;
         
        //for(k=0; k<=Ns; k++)
        //{
          //if(a_search < 0.05)
           //a_search = 0.05;
           
           // printf("a_search=%4.5e\n",a_search);
          search_resonant_frequencies(e_search, p_search, a_search , data);
          
          //a_search -= Da;
        //}
        p_search -= DP;
      }
      e_search -= DEC;
    }

  
  if( RES_COUNT == 0 )
  {
    printf("\n");
    printf("---------------------------------------------------------------------------\n\n");
    printf("[!] NO RESONANCE has been found within the selected range of parameters :( \n\n");
      
    fprintf(data,"------------------------------------------------------------------------\n\n");
    fprintf(data,"[!] No resonance has been found within the following range of parameters:\n\n");
    fprintf(data,"e_o = %4.4e, e_f = %4.4e, p_o = %4.4e, p_f = %4.4e, a = %4.4e\n\n",
                Ec_o_search, e_search, P_o_search, p_search, a_search);
    fprintf(data,"------------------------------------------------------------------------\n\n");
  }
    
 fclose(data);

    
 if(EVAL_arg == 1)
 {
    printf("\n");
    printf("In this search we are crossed invalid arguments for the Elliptical integrals\n\n");
 }

/*----- This is the end -----*/
    return 0;
}
