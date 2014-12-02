/*===============================================================================
 
 Computation of the flat cartesian coordinates, velocity and acceleration of the 
 particle
 -------------------------------------------------------------------------------
 Created by Priscilla on 19/12/2012
 -------------------------------------------------------------------------------
 Last Update : 25.09.2013
 ===============================================================================*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/types.h>
#include<time.h>
#include <math.h>

#include "global_quantities_kerr.h"
#include "global_quantities_main.h"
#include "global_quantities_osc.h"
#include "physical_quantities.h"
#include "macros_2.h"

int christoffel_symbols(double, double, double g[4][4][4]);
int dchi_dpsi_dphi_dt(double *, double, double,double *,double *,double *,double *);


int flat_space_XVA(double *x, double r, double costh, double dphidt, double drdt, double dcosthdt, double *X_bl, double *V_bl, double *A_bl)
{
  double dthdt;      // Chi_p modulus 2 pi
  double zeta;               // z = cos^2(Theta_p) = z_ cos^2(Chi_p)
  double sinth;              // cos(Theta_p), sin(Theta_p)
  double cosph,sinph;        // cos(Phi_p), sin(Phi_p)
  double r2;                 // (R_p)^2
  double drdt2, dthdt2, dphidt2;
  double D2rdt2,D2thdt2,D2phidt2;
    
    
  sinth = sqrt(1. - costh*costh);
  zeta = costh*costh;
  
  cosph = cos(x[1]);
  sinph = sqrt(1. - cosph*cosph);

  r2 = r*r;
  dthdt = -dcosthdt/sqrt(1.-zeta);
    
  drdt2 = drdt*drdt;  //check times
  dphidt2  = dphidt*dphidt;
  dthdt2   = dthdt*dthdt;
  

/*---- Cartesian coordinates from the Orbital Boyer-Lindquist coordinates----*/
  X_bl[0] = r*sinth*cosph;
  X_bl[1] = r*sinth*sinph;
  X_bl[2] = r*costh;
    
/*----- Computing the Velocities of the Euclidean Cartesian coordinates, (X,Y,Z), with respect to the Boyer-Lindquist time T_p -----*/
  V_bl[0] = drdt*sinth*cosph + r*(-costh*dcosthdt*cosph/sqrt(1.-zeta) - sinth*sinph*dphidt);
  V_bl[1] = drdt*sinth*sinph + r*(-costh*dcosthdt*sinph/sqrt(1.-zeta) + sinth*cosph*dphidt);
  V_bl[2] = drdt*costh - r*costh*dcosthdt/sqrt(1.-zeta);
    
/*----- Computing the Accelerations of the Boyer-Lindquist coordinates, (R_p,Theta_p,Phi_p), with respect to the Boyer-Lindquist time T_p -----*/
  christoffel_symbols(r,costh,GAMMA);
    
/*----- Computing the Velocities of the Boyer-Lindquist radial and polar coordinates, (Theta_p,R_p), with respect to the Boyer-Lindquist time T_p -----*/
/* 
    chi_2pi = fmod(CHI_NEW,2.0*PI);
   
   if ( GV_Flag_equatorial_orbit == 'y' )
        dthdt = 0.0;
    else
    {
        if ( (chi_2pi >= 0.0) && (chi_2pi < PI) )
        {
            dthdt = sqrt( (GV_Z_minus-zeta)/(1.0-zeta) )*GV_DChi_p_DT;
        }
        else if ( (chi_2pi >= PI) && (chi_2pi < 2.0*PI) )
        {
            dthdt = - sqrt( (GV_Z_minus-zeta)/(1.0-zeta) )*GV_DChi_p_DT;
        }
    }
    drdt = (r2/GV_P_orbit)*GV_Eccentricity*sin(GV_Psi_p_NEW)*GV_DPsi_p_DT;
*/ 
        
  D2rdt2 = GAMMA[0][1][1]*drdt2*drdt + 2.0*drdt2*GAMMA[0][1][2]*dthdt + 2.0*drdt2*GAMMA[0][1][3]*dphidt + 2.0*GAMMA[0][0][1]*drdt2
           + drdt*GAMMA[0][2][2]*dthdt2 + 2.0*drdt*dthdt*GAMMA[0][2][3]*dphidt + 2.0*drdt*GAMMA[0][0][2]*dthdt
           + drdt*GAMMA[0][3][3]*dphidt2 + 2.0*drdt*GAMMA[0][0][3]*dphidt + drdt*GAMMA[0][0][0] - GAMMA[1][1][1]*drdt2
           - 2.0*drdt*GAMMA[1][1][2]*dthdt - 2.0*drdt*GAMMA[1][1][3]*dphidt - 2.0*GAMMA[1][0][1]*drdt - GAMMA[1][2][2]*dthdt2
           - 2.0*dthdt*GAMMA[1][2][3]*dphidt - 2.0*GAMMA[1][0][2]*dthdt - GAMMA[1][3][3]*dphidt2 - 2.0*GAMMA[1][0][3]*dphidt
           - GAMMA[1][0][0];
    
  D2thdt2 = dthdt*GAMMA[0][1][1]*drdt2+2.0*drdt*GAMMA[0][1][2]*dthdt2+2.0*dthdt*drdt*GAMMA[0][1][3]*dphidt
            + 2.0*drdt*GAMMA[0][0][1]*dthdt+GAMMA[0][2][2]*dthdt2*dthdt+2.0*dthdt2*GAMMA[0][2][3]*dphidt+2.0*GAMMA[0][0][2]*dthdt2
            + dthdt*GAMMA[0][3][3]*dphidt2+2.0*dthdt*GAMMA[0][0][3]*dphidt+dthdt*GAMMA[0][0][0]-GAMMA[2][1][1]*drdt2
            - 2.0*drdt*GAMMA[2][1][2]*dthdt-2.0*drdt*GAMMA[2][1][3]*dphidt-2.0*GAMMA[2][0][1]*drdt-GAMMA[2][2][2]*dthdt2
            - 2.0*dthdt*GAMMA[2][2][3]*dphidt-2.0*GAMMA[2][0][2]*dthdt-GAMMA[2][3][3]*dphidt2-2.0*GAMMA[2][0][3]*dphidt
            - GAMMA[2][0][0];
    
  D2phidt2 = dphidt*GAMMA[0][1][1]*drdt2+2.0*dphidt*drdt*GAMMA[0][1][2]*dthdt+2.0*drdt*GAMMA[0][1][3]*dphidt2
             + 2.0*dphidt*GAMMA[0][0][1]*drdt+dphidt*GAMMA[0][2][2]*dthdt2+2.0*dthdt*GAMMA[0][2][3]*dphidt2+2.0*dphidt*GAMMA[0][0][2]*dthdt
             + GAMMA[0][3][3]*dphidt2*dphidt+2.0*GAMMA[0][0][3]*dphidt2+dphidt*GAMMA[0][0][0]-GAMMA[3][1][1]*drdt2
             - 2.0*drdt*GAMMA[3][1][2]*dthdt-2.0*drdt*GAMMA[3][1][3]*dphidt-2.0*GAMMA[3][0][1]*drdt-GAMMA[3][2][2]*dthdt2
             - 2.0*dthdt*GAMMA[3][2][3]*dphidt-2.0*GAMMA[3][0][2]*dthdt-GAMMA[3][3][3]*dphidt2-2.0*GAMMA[3][0][3]*dphidt
             - GAMMA[3][0][0];
 
/*----- Computing Accelerations of the Euclidean Cartesian coordinates, (X_p,Y_p,Z_p), with respect to the Boyer-Lindquist Time T_p -----*/
  A_bl[0] = D2rdt2*sinth*cosph + 2.0*drdt*costh*dthdt*cosph
          - 2.0*drdt*sinth*sinph*dphidt - r*sinth*dthdt*dthdt*cosph
          + r*costh*D2thdt2*cosph - 2.0*r*costh*dthdt*sinph*dphidt
          - r*sinth*cosph*dphidt*dphidt - r*sinth*sinph*D2phidt2;
    
  A_bl[1] = D2rdt2*sinth*sinph + 2.0*drdt*costh*dthdt*sinph
          + 2.0*drdt*sinth*cosph*dphidt - r*sinth*dthdt*dthdt*sinph
          + r*costh*D2thdt2*sinph + 2.0*r*costh*dthdt*cosph*dphidt
          - r*sinth*sinph*dphidt*dphidt + r*sinth*cosph*D2phidt2;
    
  A_bl[2] = D2rdt2*costh - 2.0*drdt*sinth*dthdt - r*costh*dthdt*dthdt - r*sinth*D2thdt2;
    
/*----- This is the end -----*/
    return 0;
}