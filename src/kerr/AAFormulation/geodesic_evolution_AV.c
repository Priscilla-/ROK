/*==============================================================================
 
 This routine evolves the *geodesic* equations for the Kerr Spacetime using
 angle variables wr and wth. Notice that the cooridinate time , t,is evolved 
 independently wrt the "Mino time" lambda, whereas wr and wth are evolved wrt t. 
 -------------------------------------------------------------------------------
                                            Created by Priscilla on 20/12/2012
 -------------------------------------------------------------------------------
 Last Update: 31.07.13
 =============================================================================*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"
#include "physical_quantities.h"
#include "macros_2.h"


void JacobianEllipticFunctions(double, double, double *, double * , double * );
int analytic_r_costh(double *, double *, double *);

int geodesic_evolution_AV(double *x, double *dxdt)
{
  double q3,q4,q5,q6,z;
  double sin_psi, sin_chi, cos_chi, cos_psi;
  double sigma, delta;
  double beta,omec2,omEn2;
  double pr, Tr, Tth, dtdlambda;
  double r, costh;
  
  double dchidt,numer_chi,denom_chi;
  
 /*----- mapping wr, wth into r, costh -----*/
  analytic_r_costh(x,&r, &costh);
  
 /*----- Auxiliary quantities -----*/
  omEn2 = 1.0 - SQR(ENERGY);        
  beta = SQR(SPIN)*omEn2;        
  omec2 = 1.0 - SQR(ECCENTRICITY);      
  sigma = r*r + SQR(SPIN*costh);
  delta = r*r -2*r + SPIN*SPIN;
  z = costh*costh;
  cos_psi = (P_p/r - 1.)/ECCENTRICITY;
  sin_psi = sqrt(1. - cos_psi*cos_psi);
  cos_chi = costh/sqrt(Z_minus);
    
  q3 = P_p - P3_p - ECCENTRICITY*(P_p + P3_p*cos_psi);  // p - p3 - e ( p + p3 cos(psi) )
  q4 = P_p - P4_p + ECCENTRICITY*(P_p - P4_p*cos_psi);  // p - p4 + e ( p - p4 cos(psi) )
  q5 = r*r + SQR(SPIN);     
  q6 = ENERGY*(q5*q5/delta-SQR(SPIN)) - 2.0*SPIN*r*LZ/delta;  // gamma = E [ ( r^2 + a^2 )^2 / Delta - a^2 ] - 2 M a r Lz / Delta
    
  Tr = ENERGY*( SQR(r*r + SPIN*SPIN)) - 2.0*SPIN*r*LZ;
  Tr /=delta ;

  Tth = -SPIN*SPIN*ENERGY*(1.-costh*costh);
  
/*--- Derivatives wrt observer's time ---*/
    
  dxdt[4] = Tr + Tth; // dtdlambda 
  dxdt[1] = OMEGA_r/dxdt[4];      // dwr/dt = (dwr/dlambda)*(dlambda/dt)
  
  dxdt[2] = OMEGA_th/dxdt[4];    // dwth/dt= (dwth/dlambda)*(dlambda/dt)
 
  dxdt[3] = (1.0/(1.0-z) - SQR(SPIN)/delta)*LZ + 2.0*SPIN*ENERGY*r/delta; // dphi/dt
  dxdt[3]/= q6 + SQR(SPIN)*ENERGY*z;
  
/*----- This is the end -----*/
  return 0;
}