/*==============================================================================
 This routine Computes the osculating conditions for the Orbital frequencies
 in terms of Mino time. The formulation employed here is based in the following
 references:
 (a)Gair et al [PRD 83, 044037 (2011)]: for the osculating formalism
 (b)Fujita et al [PTP, 121,4,2009]: For the action angle variable formalism
 (c)Fujita & Hikida [CQG, 26,135002]: For computing the fundamental frequencies
 ===============================================================================
 NOTE:
 The geodesic equations for bound orbits can be written employing BL coordinates
 q_µ =(t,r,theta, phi). And for being a bound orbit, the geodesics have
 associated tree fundamental frequencies Omega_r, Omega_thta, Omeg_phi, in terms
 of the observer's time. This frequencies give the multiply-periodic motion in r,
 theta and phi. However,the motion in r and theta is not quite periodic because
 both movements are coupled, which makes hard to find the coefficients of any
 fourier expansion based on these fundmantal frequencies.
 However, in tems of the "Mino time", lambda, the motion is decoupled and we can
 find angle variables w_r = omega_r*t=Omega_r*λ, w_theta = omega_theta*t
 =Omega_theta*λ, such that:
 r(Ω_r) and theta(Ω_theta), and computing the Furier coefficients is
 easier.
 
 On the other hand, the osculatingevolution conditions:
 
                (∂q_µ/∂I^µ)*(dI^µ/dtau)=0,
 where
            I^µ = {I^A,I^B},  with I^A={E,L_z,Q} and I^B={W^r,w^th},
 
 are the orbital parameters, leads to the following equations:
 
 ∂r(Ω_r)/∂I^µ *dI^µ/dtau= ∂r/∂Ω_r*∂Ω_r/∂I^µ *dI^µ/dtau = 0,
 
 so eventually we just have to solve for:
 
        ∂I^B/dtau = -(∂q_µ/∂I^b)^-1(∂q_µ/∂I^A)(dI^A/dtau)
 
 -------------------------------------------------------------------------------
 Created by Priscilla on 26/03/2013
 -------------------------------------------------------------------------------
 Last Update: 08/04/2013                    Modified by PCM
 =============================================================================*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"

#include "macros_1.h"
#include "macros_2.h"

double ellK(double, double);

void JacobianEllipticFunctions(double, double, double *, double * , double * );
int osculating_evolution(double, double, double, double, double, double, double, double, double, double, double, double,double*,double*);
int PN_ELzQdot( double *, double *, double *);

int RHs_r_theta_phi(double *drdl, double *dcosthdl, double *dphidl, double rp, double cosTh)
{
  double sn_r, sn_th;                     // Jacobi's elliptic functions
  double cn_r, cn_th;                     // Jacobi's elliptic functions
  double dn_r, dn_th;                     // Jacobi's elliptic functions
  double varphi_r,varphi_th;              // Argument Jacobi's elliptic functions
  double drdw, dcosthdw;                  // Radial and plolar derivatives wrt their fundamental frequencies
  double dwrdL, dwthdL;                   // Osculating evolution for the fundamental frequencies
  double denom;                           // Aux variables
  double r,r2pa2,delta,SQr2pa2,gamma,z;
  double denom_phi,numer_phi;
  double cos_psi;
  double W_r, W_th;                       // Action Angle variables
  int sign_r, sign_th;                    // Aux variables

    double GV_E = GV_E;
    double GV_Lz = GV_Lz;
    
    
/*---Computing fundamental frequencies for a distant observer---*/
    w_obs_r = GV_Omega_r/GV_Omega_t;
    w_obs_th = GV_Omega_th/GV_Omega_t;
    
/*---Action Angle variables ----*/
    W_r  = w_obs_r*dt;
    W_th = w_obs_th*dt;

/*----- Auxiliary quantities -----*/
    r2pa2= SQR(rp) + GV_Spin2;                        // r^2 + a^2
    SQr2pa2 = SQR(r2pa2);                             // (r^2+a^2)^2
    delta =  r2pa2 -2*rp;
    gamma = GV_E*(SQr2pa2/delta-GV_Spin2) - 2.0*GV_Spin*rp*GV_Lz/delta; // gamma = E [ ( r^2 + a^2 )^2 / Delta - a^2 ] - 2 M a r Lz / Delta
    
    numer_phi = ( (1.0/(1.0-SQR(cosTh)) - GV_Spin2/delta)*GV_Lz + 2.0*GV_Spin*GV_E*rp/delta );
    denom_phi = sqrt((GV_Spin2*(1.0 - SQR(GV_E)))*(GV_Z_plus-SQR(cosTh)));
    
/*---Arguments of the Jacobi elliptic functions----*/
    if( (W_r>= 0)&&(W_r<=PI) )
    {
        varphi_r = W_r*GV_Kcomp_r/PI;
        sign_r = 1;
    }
    else if( (W_r>= PI)&&(W_r<=2*PI) )
    {
        varphi_r = (2*PI-W_r)*GV_Kcomp_r/PI;
        sign_r = -1;
    }
    
    if( (W_th>= 0)&&(W_th<=PI/2) )
    {
        varphi_th = W_th*2*GV_Kcomp_th/PI;
        sign_th = 1;
    }
    else if( (W_th>= PI/2)&&(W_th<=3*PI/2) )
    {
        varphi_th = (PI-W_th)*2*GV_Kcomp_th/PI;
        sign_th =-1;
    }
    else if( (W_th>= 3*PI/2)&&(W_th<=PI) )
    {
        varphi_th = (W_th-2*PI)*2*GV_Kcomp_th/PI;
        sign_th = 1;
    }
    
/*---Jacobi elliptic functions----*/
    JacobianEllipticFunctions(varphi_r, SQR(GV_kr),  &sn_r,  &cn_r ,  &dn_r);
    JacobianEllipticFunctions(varphi_th,SQR(GV_kth), &sn_th, &cn_th , &dn_th);
    
/*--- Radial derivative wrt Wr ----*/
    denom = SQR(sn_r)*(GV_R_apo - GV_R_peri)-(GV_R_apo - GV_R_3);
    
    drdw = 2*sign_r*sn_r*cn_r*dn_r*(GV_R_apo - GV_R_peri)*(GV_R_apo - GV_R_3)*(GV_R_peri - GV_R_3)*GV_Kcomp_r;
    drdw /= SQR(denom)*PI;

/*--- Polar derivative wrt Wth ----*/
    dcosthdw = 2*sign_th*sqrt(GV_Z_minus)*cn_th*dn_th*GV_Kcomp_th;
    dcosthdw /= PI;

/*--- Evolution along a geodesic ---*/
 
    //--Radial geodesic equation 
     *drdl = drdw*GV_Omega_r;
      
    //--Polar geodesic equation
     *dcosthdl = dcosthdw*GV_Omega_th;
  
    //--Azimutal geodesic equation
    *dphidl = numer_phi/denom_phi;

/*----- This is the end -----*/
    return 0;
}