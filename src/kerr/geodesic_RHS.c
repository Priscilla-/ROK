/*==============================================================================
   This Routine evaluates the Right-Hand Sides (RHSs) of the system of ODEs that
   governs the dynamics of geodesics of a Kerr Black Hole in Boyer-Lindquist
   coordinates: dPsi/dt, dChi/dt, dPhi/dt.[dXdt, with t BL time]
   In addition it evaluates the RHS of the "Mino" time Lambda: dlambda/dt
-------------------------------------------------------------------------------
                                             Created by Priscilla on 19/12/2012
-------------------------------------------------------------------------------
Last Update : 24.07.13
==============================================================================*/

#include<stdio.h>
#include<stdlib.h>

#include "global_quantities_kerr.h"
#include "global_quantities_main.h"
#include "macros_2.h"


int geodesic_RHS(double *z_chi, double *z_psi, double *z_phi,double *z_lambda, double chi, double psi)
{
  double q1,q2,q3,q4,q5,q6,delta,z1;
  double denom_psi,denom_phi,denom_chi,numer_psi,numer_phi,numer_chi;
  double cos_psi,sin_psi,cos_chi;
  double beta,omec2,omEn2;
  
/*----- Auxiliary quantities -----*/
  omEn2 = 1.0 - SQR(ENERGY);      // 1 - E^2
  beta = SQR(SPIN)*omEn2;         // a^2 (1 - E^2)
  omec2 = 1.0 - SQR(ECCENTRICITY);// 1 - e^2
    
  cos_psi = cos(psi);   
  sin_psi = sin(psi);        
  cos_chi = cos(chi);
 
  
  q1 = 1.0 + ECCENTRICITY*cos_psi;
  q2 = P_p/q1;                     // r = p M / ( 1 + e cos(psi) )
  q3 = P_p - P3_p - ECCENTRICITY*(P_p + P3_p*cos_psi);  // p - p3 - e ( p + p3 cos(psi) )
  q4 = P_p - P4_p + ECCENTRICITY*(P_p - P4_p*cos_psi);  // p - p4 + e ( p - p4 cos(psi) )

  q5 = SQR(q2) + SQR(SPIN);   // r^2 + a^2
  delta = q5 - 2.0*q2;        // Delta = r^2 - 2 M r + a^2
 
  q6 = ENERGY*(SQR(q5)/delta-SQR(SPIN)) - 2.0*SPIN*q2*LZ/delta;  // gamma = E [ ( r^2 + a^2 )^2 / Delta - a^2 ] - 2 M a r Lz / Delta
						               
  z1 = Z_minus*SQR(cos_chi);            // z = z_minus cos(chi)^2 = cos(theta)^2

  numer_chi = sqrt(beta*(Z_plus-z1));   // dchi/dlambda = sqrt( Beta*(Z_plus - z) )
  denom_chi = q6 + SQR(SPIN)*ENERGY*z1; // dt/dlambda
  
  numer_psi = sqrt(omEn2*q3*q4);
  denom_psi = omec2*denom_chi;
  
  numer_phi = ( (1.0/(1.0-z1) - SQR(SPIN)/delta)*LZ + 2.0*SPIN*ENERGY*q2/delta );
  denom_phi = denom_chi;
  
/*----- RHS Evolution equation for Chi -----*/
  if ( FLAG_equatorial_orbit == 'y' )
    *z_chi = 0.0;
  else
    *z_chi = numer_chi/denom_chi;

/*----- RHS Evolution equation for Psi -----*/
  if ( FLAG_circular_orbit == 'y' )
    *z_psi = 0.0;
  else
    *z_psi = numer_psi/denom_psi;
  
/*----- RHS Evolution equation for Phi -----*/
  *z_phi = numer_phi/denom_phi;

/*----- RHS Evolution equation for Lamda ("Mino time") -----*/
  *z_lambda = 1./denom_chi;
  
  
/*-----This is the end -----*/
  return 0;

}

