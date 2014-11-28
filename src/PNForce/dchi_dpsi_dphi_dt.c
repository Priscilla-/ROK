/*==============================================================================
   This Routine evaluates the Right-Hand Sides (RHSs) of the system of ODEs that
   governs the dynamics of geodesics of a Kerr Black Hole in Boyer-Lindquist
   coordinates.
-------------------------------------------------------------------------------
                                             Created by Priscilla on 19/12/2012
-------------------------------------------------------------------------------
Last Update : 10.12.12
==============================================================================*/

#include<stdio.h>
#include<stdlib.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"
#include "macros_2.h"
#include "physical_quantities.h"

int analytic_r_costh(double *, double *, double *);

int dchi_dpsi_dphi_dt(double *dxdt, double *x, double r, double costh, double *dpsidt, double *dchidt,double *dphidt, double *dtaudt)
{
  double r2pa2,SQr2pa2,gamma,z;
  double q3,q4,q5;
  double dpsidlambda, dchidlambda, dphidlambda, dtdlambda;
  double sin_psi, sin_chi, cos_chi, cos_psi;
  double sigma, delta;
  double drdtau, dcosthdtau, dzmdt;
  double beta,omec2,omEn2;

/*----- Auxiliary quantities -----*/
  cos_psi = (x[3]/r - 1.)/x[4];
  cos_chi = costh/sqrt(Z_minus);

  sin_psi = sqrt(1. - cos_psi*cos_psi);

  omEn2 = 1.0 - SQR(ENERGY);  // 1 - E^2
                              //beta = SQR(SPIN)*omEn2;     // a^2 (1 - E^2)
  omec2 = 1.0 - x[4]*x[4];    // 1 - e^2
  
  
  q3 = x[3] - P3_p - x[4]*(x[3] + P3_p*cos_psi);  // p - p3 - e ( p + p3 cos(psi) )
  q4 = x[3] - P4_p + x[4]*(x[3] - P4_p*cos_psi);  // p - p4 + e ( p - p4 cos(psi) )
  q5 = r*r + SQR(SPIN);
  
  sigma = r*r + SPIN*SPIN*costh*costh;
  delta = q5 - 2.0*r;        // Delta = r^2 + a^2 - 2 M r
  
  dpsidlambda = sqrt(omEn2*q3*q4)/omec2;
  dchidlambda = sqrt(SPIN*SPIN*omEn2*(Z_plus -costh*costh));
  dphidlambda = (ENERGY*q5 -SPIN*LZ)*SPIN/delta + LZ/(1.-costh*costh)-SPIN*ENERGY;
  dtdlambda   = (ENERGY*q5 -SPIN*LZ)*q5/delta -SPIN*SPIN*ENERGY*(1.-costh*costh) + SPIN*LZ;
  

/*----- RHS Evolution equation for Psi -----*/
  if ( FLAG_circular_orbit == 'y' )
   *dpsidt = 0.0;
  else
   *dpsidt = dpsidlambda/dtdlambda;

/*----- RHS Evolution equation for Chi -----*/
  if ( FLAG_equatorial_orbit == 'y' )
    *dchidt = 0.0;
  else
    *dchidt =dchidlambda/dtdlambda;

/*----- RHS Evolution equation for Phi -----*/
  *dphidt = dphidlambda/dtdlambda;

/*----- dtaudt -----*/
  *dtaudt = sigma/dtdlambda;//( (ENERGY*(q5*q5/delta - SQR(SPIN)) - 2.*q2*SPIN*LZ/delta + SQR(SPIN)*ENERGY*z1));

/*-----This is the end -----*/
  return 0;

}

