/*==============================================================================
 
            This file reads the physical and computational parameters

--------------------------------------------------------------------------------
                                            Created by Priscilla on 19/12/2012
--------------------------------------------------------------------------------
 Last Update : 11.12.12
==============================================================================*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"

int reading_parameters(void)
{
  FILE *data;
 // int i;
  char buffer[10000];

/*----- Reading Evolution Parameters -----*/
  data=fopen("parameters","r");
      fscanf(data, "%[^\n]\n", buffer);
      fscanf(data, "%[^\n]\n", buffer);
      fscanf(data, "%[^\n]\n", buffer); 
    fscanf(data, "%*90c %lf\n", &EMRI_Parameter[1]);     // Black Hole Mass: M. (M_sun units)
    fscanf(data, "%*90c %lf\n", &EMRI_Parameter[2]);     // Black Hole Spin: S.  (M. units)
    fscanf(data, "%*90c %lf\n", &EMRI_Parameter[3]);     // Mass Ratio: mu
    fscanf(data, "%*90c %lf\n", &EMRI_Parameter[4]);     // Eccentricity of the Orbit at the Initial Time, t = t_o (M. units): e_o
    fscanf(data, "%*90c %lf\n", &EMRI_Parameter[5]);     // Dimensionless Semi-Latus Rectum of the Orbit at the Initial Time, t = t_o (M. units): p_o
    fscanf(data, "%*90c %lf\n", &EMRI_Parameter[6]);     // Inclination of the Orbit at the Initial Time, t = t_o (M. units): theta_inc_o
    fscanf(data, "%*90c %lf\n", &EMRI_Parameter[7]);     // EMRI Polar Angle with respect to the Ecliptic Reference Frame: theta_S
	fscanf(data, "%*90c %lf\n", &EMRI_Parameter[8]);       // EMRI Azimuthal Angle with respect to the Ecliptic Reference Frame: phi_S
    fscanf(data, "%*90c %lf\n", &EMRI_Parameter[9]);     // Spin Vector Polar Angle with respect to the Ecliptic Reference Frame: theta_K
	fscanf(data, "%*90c %lf\n", &EMRI_Parameter[10]);      // Spin Vector Azimuthal Angle with respect to the Ecliptic Reference Frame: phi_K
    fscanf(data, "%*90c %lf\n", &EMRI_Parameter[11]);    // Distance to the EMRI from the Ecliptic Baricentric Reference Frame D_L divided by the mass ratio mu (Gpc units)
    fscanf(data, "%*90c %lf\n", &EMRI_Parameter[12]);    // Boyer-Lindquist Angle psi [associated with r_BL] at t = t_o: psi_o
	fscanf(data, "%*90c %lf\n", &EMRI_Parameter[13]);    // Boyer-Lindquist Angle chi [associated with theta_BL] at t = t_o: chi_o
	fscanf(data, "%*90c %lf\n", &EMRI_Parameter[14]);    // Boyer-Lindquist Angle varphi [identical to varphi_BL] at t = t_o: varphi_o
       fscanf(data, "%[^\n]\n", buffer);
       fscanf(data, "%[^\n]\n", buffer);
       fscanf(data, "%[^\n]\n", buffer);
    fscanf(data, "%*82c %c\n", &FLAG_resonances);        // Flag to control whether we are searching for resonant orbits
    fscanf(data, "%*82c %u %u \n",&Resorder[0],&Resorder[1]); // Resonance order (k*w_r - n*w_th=0): Resorder[0] = k, Resorder[1] = n
    fscanf(data, "%*82c %lf\n", &Res_accuracy);          // Resonance condition accuracy (k*w_th - n*w_r)
    fscanf(data, "%*82c %lf\n", &EMRI_Parameter[15]);    // Final eccentricity in resonance searches
    fscanf(data, "%*82c %lf\n", &Real_Parameter[6]);     // Eccentricity step-size in resonance searches
    fscanf(data, "%*82c %lf\n", &EMRI_Parameter[16]);    // Final semilatus-rectum in resonance seraches
    fscanf(data, "%*82c %lf\n", &Real_Parameter[7]);     // Semilatus-rectum step-size in resonance searches
      fscanf(data, "%[^\n]\n", buffer);
      fscanf(data, "%[^\n]\n", buffer);
      fscanf(data, "%[^\n]\n", buffer);
    fscanf(data, "%*82c %c %lf\n", &FLAG_OS, &MAG_SF);    // Use osculating evolution? include post-Newtonian SF? 
    fscanf(data, "%*82c %c\n", &FLAG_incConsPiece);       // Include conservative Effects? [y/n]
    fscanf(data, "%*82c %c\n", &FLAG_RR );                // Apply Radiation Reaction Effects? [y/n]
     fscanf(data, "%*82c %lf\n",&STORE_rate);             // Flag to store the EMRI evolution
      fscanf(data, "%[^\n]\n", buffer);
      fscanf(data, "%[^\n]\n", buffer);
    fscanf(data, "%*82c %lf\n", &TIME_yrs);                  // Evolution Time [yrs]
	fscanf(data, "%*82c %lf\n", &DTsec);                 // Sampling Time/Data Saving Time Interval [sec]
	fscanf(data, "%*82c %lf\n", &DT_INI);                // ODE Solver Starting Time Step [Black Hole Mass]
	fscanf(data, "%*82c %lf\n", &DT_MIN);                // ODE Solver Minimum Time Step [Black Hole Mass]
	fscanf(data, "%*82c %lf\n", &ACCURACY);              // ODE Solver Accuracy [Black Hole Mass]
    fscanf(data, "%[^\n]\n", buffer);
    fscanf(data, "%[^\n]\n", buffer);
  fclose(data);
   
// printf(" %4.6e\n",STORE_rate);exit(0);
    
 /*---- Warnings on the screen ---*/
    if( Resorder[1] > Resorder[0] )
    {
        printf(" [!] Resonance order: k should be smaller than or equal to n. Correct this and have a coffee! \n\n");
        exit(0); 
    }
    
    if(( EMRI_Parameter[4] < EMRI_Parameter[15] )||( EMRI_Parameter[5] < EMRI_Parameter[16] ))
    {
        printf(" [!] The convention in ROK is e_o > e_f and p_o > p_f. Change it accordingly \n\n");
        exit(0);
    }
    
    if((FLAG_OS =='y')&& (FLAG_RR=='y') )
      printf(" [!] You have chosen two types of evolution: one with PN sefl-force and one with Radiative effects! \n\n");
        

/*----- This is the end -----*/
  return 0;

}

