/*==============================================================================
 
 This routine computes the fundamental frequencies associated with a given
 geodesic in "Mino" time λ: Omega_r, Omega_theta, Omega_phi, Omeag_t.
 Where the angular frequencies of the radial and the polar motion are:
 
 Omega_r = 2PI/P_r
 
 Omega_theta =2PI/P_theta
 
 with P_r and P_theta the fundamental periods for the radial and the polar
 motions with respect to λ.
 
 From these expressions we can also obtain the frequencies in terms of an
 observer's time t:
 
 omega_r = Omega_r/Omeag_t
 
 omega_phi = Omega_phi/Omeag_t
 
 omega_th = Omega_th/Omega_t;
 
 
 The expressions computed here are quoted in Fujita & Hikida
 (Class. Quantum Grav. 26 (2009))
 -------------------------------------------------------------------------------
 Created by Priscilla on 20/12/2012
 -------------------------------------------------------------------------------
 Last Update : 26.04.13
 ==============================================================================*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "physical_quantities.h"
#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"

#include "macros_2.h"

double ellK(double, double);
double ellE(double, double);
double ellP(double, double, double );

int compute_frequencies(double *FM_W_r, double *FM_W_th, double *FM_W_t, double *K_r, double *K_th,
                        double *E_r, double *E_th, double *S_kr, double *S_kth)
{
  double Kcomp_r, Kcomp_th;                                     // Complete elliptic integral of the firs kind
  double Pcomp_mhp, Pcomp_phm, Pcomp_mhr, Pcomp_mhm, Pcomp_mZm; // Complete elliptic integral of the third kind
  double Ecomp_r, Ecomp_th;                                     // Complete elliptic integral of the second kind
  double kr, kr_num, kr_dem;                                    // Arguments for the complete elliptic integrals
  double kth;                                                   // Arguments for the complete elliptic integrals
  double hr, hp, hm;                                            // Arguments for the complete elliptic integrals
  double rp, rm;                                                // Auxiliar variables: roots of Delta = r^2-2*Mr+a^2 NOTE: M = 1 here!!
  double q1p, q1m, q2p, q2m, g1, g2, g3, g4p, g4m, g5, g6, g7;  // Auxiliar variables used in computing Omgega_phi
  double M_Omega_th, M_Omega_r, M_Omega_phi, M_Omega_t;         // Fundamental frequencies associated with Mino time
  double omega_r, omega_th, omega_phi;                          // Fundamental frequencies seen from a distant observer
  double mode_th, mode_r;                                       // Accuracies
  
   
/*---Computing the arguments of the elliptical integrals---*/
  kr_num = (R_apo - R_peri)*(R_3 - R_4);
  kr_dem = (R_apo - R_3)*(R_peri - R_4);
  mode_r = kr_num/kr_dem;
  kr = sqrt(mode_r);

  mode_th = Z_minus/Z_plus;
  kth = sqrt(mode_th);
  
  if((R_apo < R_peri)||(R_3 < R_4))
  {  printf("Swapped turning points\n");
    exit(0);
  }
  
  rp = 1.0 + sqrt(1.0-SQR(SPIN)); // rp = M + [M^2 - a^2]^(1/2)
  rm = 1.0 - sqrt(1.0-SQR(SPIN)); // rm = M - [M^2 - a^2]^(1/2)
    
  hr = (R_apo-R_peri)/(R_apo-R_3);
  hp = hr*(R_3 - rp)/(R_peri - rp);
  hm = hr*(R_3 - rm)/(R_peri - rm);
  

/*---Computing Legendre elliptical integrals First kind---*/
  Kcomp_r  = ellK(0.5*PI,kr);
  Kcomp_th = ellK(0.5*PI,kth);
    
/*---Computing Legendre elliptical integrals Second kind---*/
  Ecomp_r =  ellE(0.5*PI,kr);
  Ecomp_th = ellE(0.5*PI,kth);
    
/*---Computing Legendre elliptical integrals Third kind---*/    
  Pcomp_mhp = ellP(0.5*PI, -hp, kr);
  Pcomp_mhm = ellP(0.5*PI, -hm, kr);
  Pcomp_phm = ellP(0.5*PI,  hm, kr);
  Pcomp_mhr = ellP(0.5*PI, -hr, kr);
  Pcomp_mZm = ellP(0.5*PI, -Z_minus, kth);
    
/*---Computing fundamental frequencies in 'Mino time'---*/
/**---auxiliar variables----**/
  q1p = (2.0*ENERGY*rp - SPIN*LZ)/(R_3 - rp);
  q1m = (2.0*ENERGY*rm - SPIN*LZ)/(R_3 - rm);
    
  q2p = Kcomp_r - (R_peri - R_3)/(R_peri - rp)*Pcomp_mhp;
  q2m = Kcomp_r - (R_peri - R_3)/(R_peri - rm)*Pcomp_mhm;
    
  g1 = 2.0*SPIN*ENERGY*sqrt(Z_plus)*(Kcomp_th - Ecomp_th)/(PI*sqrt(1.0 - SQR(ENERGY) ));
    
  g2 = 0.5*ENERGY*((R_3*(R_apo + R_peri + R_3) - R_apo*R_peri)*Kcomp_r
     + (R_peri - R_3)*(R_apo + R_peri + R_3 +R_4)*Pcomp_mhr
     + (R_apo-R_3)*(R_peri-R_4)*Ecomp_r);
    
  g3 = 2.0*ENERGY*(R_3*Kcomp_r + (R_peri-R_3)*Pcomp_mhr);
    
  g4p = 2.0*( ((4.0*ENERGY - SPIN*LZ)*rp - 2.0*SQR(SPIN)*ENERGY)*(Kcomp_r - (R_peri-R_3)*Pcomp_mhp/(R_peri - rp) )/(R_3 - rp) )/(rp-rm);
  g4m = 2.0*( ((4.0*ENERGY - SPIN*LZ)*rm - 2.0*SQR(SPIN)*ENERGY)*(Kcomp_r - (R_peri-R_3)*Pcomp_mhm/(R_peri - rm) )/(R_3 - rm) )/(rp-rm);
  
  g5 = 2.0/(PI*sqrt((1.0 - SQR(ENERGY))*(R_apo - R_3)*(R_peri - R_4)));
  g6 = 2.0*SPIN/(PI*(rp - rm)*sqrt((1.0 - SQR(ENERGY) )*(R_apo - R_3)*(R_peri - R_4)) );
  g7 = 2.0*LZ/(PI*SPIN*sqrt((1.0 - SQR(ENERGY))*Z_plus));
    
/**---- Fundamental frequencies associated with Mino time----**/
  M_Omega_r  = 0.5*PI*sqrt((1.0 - SQR(ENERGY) )*(R_apo - R_3)*(R_peri - R_4))/(Kcomp_r);

  M_Omega_th = 0.5*PI*SPIN*sqrt((1.0 - SQR(ENERGY))*Z_plus)/(Kcomp_th);
    
  M_Omega_phi = g7*Pcomp_mZm*M_Omega_th + ( q1p*q2p - q1m*q2m)*g6*M_Omega_r ;
    
  M_Omega_t = 4.0*ENERGY + g1*M_Omega_th + ( g2  + g3 + g4p - g4m)*g5*M_Omega_r;
    
/*---Fundamental frequencies (in coordinate time) associated to a distant observer's time---*/
  omega_r   = M_Omega_r/M_Omega_t;
  omega_phi = M_Omega_phi/M_Omega_t;
  omega_th  = M_Omega_th/M_Omega_t;
    
/*---- Variables needed in other parts of the code ---*/
  *FM_W_r   = M_Omega_r;
  *FM_W_th  = M_Omega_th;
  *FM_W_t   = M_Omega_t;
    
  *K_r  = Kcomp_r;
  *K_th = Kcomp_th;
    
  *E_r  = Ecomp_r;
  *E_th = Ecomp_th;
    
  *S_kr = kr;
  *S_kth = kth;
    
    // Sanity check: comparing with Fujita's frequencies
    // a = 0.9
   //  printf(" omega_r = %4.6e, omega_th = %4.6e, omega_phi = %4.6e, error_r= %4.6e , error_th= %4.6e , error_phi= %4.6e \n\n", omega_r, omega_th, omega_phi, (omega_r-1.8928532285101992e-2)/(1.8928532285101992e-2), (omega_th-2.7299110395017517e-2)/(2.7299110395017517e-2), (omega_phi-3.0550463796964692e-2)/(3.0550463796964692e-2));
    // a = 1 Need to be implemented :(
  //  printf(" omega_r = %4.6e, omega_th = %4.6e, omega_phi = %4.6e, error_r= %4.6e , error_th= %4.6e , error_phi= %4.6e \n\n", omega_r, omega_th, omega_phi, (omega_r-1.9343466898960462e-2)/(1.9343466898960462e-2), (omega_th-2.6337035996626332e-2)/(2.6337035996626332e-2), (omega_phi-2.9662029663040452e-2)/(2.9662029663040452e-2));

 
       
/*----- This is the end -----*/
  return 0;
}
        

























