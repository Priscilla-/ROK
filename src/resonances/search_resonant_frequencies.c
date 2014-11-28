/*==============================================================================
 
 This routine searches the fundamental frequencies associated with a 
 given resonant geodesic (with a fixed k*w_th + n*w_r condition), once the 
 orbital parameters (e,p,inc) are specified.
 The expressions computed here are quoted in Fujita & Hikida
 (Class. Quantum Grav. 26 (2009))
 -------------------------------------------------------------------------------
                                            Created by Priscilla on 30/12/2012
 -------------------------------------------------------------------------------
 Last Update : 19.01.13  by PCM
 =============================================================================*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"

#include "macros_2.h"

double ellK(double, double);
double ellE(double, double);
double ellP(double, double, double );

int pez_to_ELCQ(char, double, double, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int turning_points(double, double, double, double, double, double, double *, double *, double *, double *, double *);

int search_resonant_frequencies(double e_res, double p_res, double a_res, FILE *data)
{
  double Kcomp_r, Kcomp_th;                                     // Complete elliptic integral of the firs kind
  double Pcomp_mhp, Pcomp_phm, Pcomp_mhr, Pcomp_mhm, Pcomp_mZm; // Complete elliptic integral of the third kind
  double Ecomp_r, Ecomp_th;                                     // Complete elliptic integral of the second kind
  double kr, kr_num, kr_dem;                                    // Arguments for the complete elliptic integrals
  double kth;                                                   // Arguments for the complete elliptic integrals
  double hr, hp, hm;                                            // Arguments for the complete elliptic integrals
  double rp, rm;                                                // Auxiliar variables: roots of Delta = r^2-2*Mr+a^2 NOTE: M = 1 here!!
  double q1p, q1m, q2p, q2m, g1, g2, g3, g4p, g4m, g5, g6, g7;  // Auxiliar variables used in computing Omgega_phi
  double mode_th, mode_r;                                       // Accuracies
  double inclination, r_peri, r_apo;
  double theta_minus, z_minus;
  double Energy, L_z, C_carter, Q_carter;
  double r_3, r_4, p_3, p_4, z_plus;
  double Omega_r, Omega_th, Omega_t;
  double resonance_condition;          // k*w_th + n*w_r = 0
    
/*-----Computing  New Orbital Parameters -----*/
    pez_to_ELCQ(FLAG_circular_orbit,  p_res, e_res, THETA_inc, a_res, &inclination, &r_peri,
                &r_apo, &theta_minus, &z_minus,&Energy, &L_z, &C_carter, &Q_carter);
    
/*----- Checking whether the Carter constant C has become negative -----*/
    if ( Q_carter < 0.0 ){printf("\n"); printf("[!]We've reached Q < 0\n"); printf("\n"); exit(0);}
    
/*-----Computing the turning points of the Boyer-Linquist radial coordinate 'r' and 'theta' -----*/
    turning_points( p_res, e_res, Energy, L_z, Q_carter, z_minus, &r_3, &r_4, &p_3, &p_4, &z_plus);
 
    
/*---Computing the arguments of the elliptical integrals---*/
    kr_num = (r_apo - r_peri)*(r_3 - r_4);
    kr_dem = (r_apo - r_3)*(r_peri - r_4);
    mode_r = kr_num/kr_dem;
    kr = sqrt(mode_r);
    
    mode_th = z_minus/z_plus;
    kth = sqrt(mode_th);
    
    rp = 1.0 + sqrt(1.0-SQR(SPIN)); // rp = M + [M^2 - a^2]^(1/2)
    rm = 1.0 - sqrt(1.0-SQR(SPIN)); // rm = M - [M^2 - a^2]^(1/2)
    
    hr = (r_apo-r_peri)/(r_apo-r_3);
    hp = hr*(r_3 - rp)/(r_peri - rp);
    hm = hr*(r_3 - rm)/(r_peri - rm);
    
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
    Pcomp_mZm = ellP(0.5*PI, -z_minus, kth);
    
/*---Computing fundamental frequencies in 'Mino time'---*/
/**---auxiliar variables----**/
    q1p = (2.0*Energy*rp - SPIN*L_z)/(r_3 - rp);
    q1m = (2.0*Energy*rm - SPIN*L_z)/(r_3 - rm);
    
    q2p = Kcomp_r - (r_peri - r_3)/(r_peri - rp)*Pcomp_mhp;
    q2m = Kcomp_r - (r_peri - r_3)/(r_peri - rm)*Pcomp_mhm;
    
    g1 = 2.0*SPIN*Energy*sqrt(z_plus)*(Kcomp_th - Ecomp_th)/(PI*sqrt(1.0 - SQR(Energy) ));
    
    g2 = 0.5*Energy*((r_3*(r_apo + r_peri + r_3) - r_apo*r_peri)*Kcomp_r
                     + (r_peri - r_3)*(r_apo + r_peri + r_3 +r_4)*Pcomp_mhr
                     + (r_apo-r_3)*(r_peri-r_4)*Ecomp_r);
    
    g3 = 2.0*Energy*(r_3*Kcomp_r + (r_peri-r_3)*Pcomp_mhr);
    
    g4p = 2.0*( ((4.0*Energy - SPIN*L_z)*rp - 2.0*SQR(SPIN)*Energy)*(Kcomp_r - (r_peri-r_3)*Pcomp_mhp/(r_peri - rp) )/(r_3 - rp) )/(rp-rm);
    g4m = 2.0*( ((4.0*Energy - SPIN*L_z)*rm - 2.0*SQR(SPIN)*Energy)*(Kcomp_r - (r_peri-r_3)*Pcomp_mhm/(r_peri - rm) )/(r_3 - rm) )/(rp-rm);
    
    g5 = 2.0/(PI*sqrt((1.0 - SQR(Energy))*(r_apo - r_3)*(r_peri - r_4)));
    g6 = 2.0*SPIN/(PI*(rp - rm)*sqrt((1.0 - SQR(Energy) )*(r_apo - r_3)*(r_peri - r_4)) );
    g7 = 2.0*L_z/(PI*SPIN*sqrt((1.0 - SQR(Energy))*z_plus));
    
/**----Radial and polar Fundamental frequencies ----**/
    Omega_r  = 0.5* PI*sqrt((1.0 - SQR(Energy) )*(r_apo - r_3)*(r_peri - r_4))/(Kcomp_r);
    
    if(Omega_r == 0.0)
    {
      printf("[!] We have hit the separatrix");
      exit(0);
    }

    Omega_th = 0.5*PI*SPIN*sqrt((1.0 - SQR(Energy))*z_plus)/(Kcomp_th);
    Omega_t  = 4.0*ENERGY + g1*Omega_th + ( g2  + g3 + g4p - g4m)*g5*Omega_r;
    
/**----Resonant condition (notice that it is given in terms of the orbital frequencies) ----**/    
    resonance_condition = fabs( Resorder[1]*Omega_th - Resorder[0]*Omega_r );
    
     
    if (resonance_condition <= Res_accuracy)
    {
        printf("---------------------------------------------------------------------------\n\n");
        printf("[*] RESONANCE! [n = %d, k = %d] n*W_th-k*W_r = %4.12e \n\n", Resorder[1], Resorder[0],resonance_condition );
        printf("[*] Frequencies in 'Mino' time: \n");
        printf(" 	 Omega_th= %4.12e  Omega_r= %4.12e \n", Omega_th,  Omega_r);
        printf("[*] Orbital parameters:\n");
        printf("  e= %4.6e   p = %4.6e  Spin = %4.6e \n\n", e_res,  p_res, a_res);

        RES_COUNT = 1;
        
        OMEGA_r_res = Omega_r;
        OMEGA_th_res = Omega_th;
        
     
        fprintf(data,"%d   %d  %4.6e  %4.10e   %4.12e   %4.12e   %4.12e %4.12e \n\n",
                Resorder[0], Resorder[1], OMEGA_r_res, OMEGA_th_res,e_res,  p_res, a_res, resonance_condition);
    
     }
/*----- This is the end -----*/
    return 0;
}

