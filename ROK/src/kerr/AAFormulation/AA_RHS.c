/*==============================================================================
 This routine Computes the osculating conditions for the Orbital frequencies
 defined terms of Mino time. The formulation employed here is based in the
 following references:
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
 find angle variables Wr = omega_r*t=Omega_r*λ, Wtheta = omega_theta*t
 =Omega_theta*λ, such that:
 r(Ω_r) and theta(Ω_theta), and computing the Furier coefficients is
 easier.
 
 On the other hand, the osculating evolution conditions given by:
 
 (∂q_µ/∂I^µ)*(dI^µ/dtau)=0,
 
 where, abusing of notation, q_µ = (r, cos(theta)) and
 
 I^µ = {I^A,I^B},  with I^A={E,L_z,Q} and I^B={W^r_o,w^th_o},
 
 lead to the following equations:
 
 ∂I^B/dλ = -(∂q_µ/∂I^A)(dI^A/dλ)/[∂q_µ/∂I^b]
 
 which diverges at ∂q_µ/∂I^b= 0. In order to avoid these turning points we also
 make use of the following osculation equation:
 
 ∂I^B/∂λ = ∂q_µ/∂λ*[F_µ +2*∂(∂q_µ/∂λ)/∂I^a]dI^a/dλ]/(∂V_µ/∂q_µ)
 
 -------------------------------------------------------------------------------
                                                        Priscilla on 12/05/2013
 -------------------------------------------------------------------------------
 Last Update: 29/09/2013                                      
 =============================================================================*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h> 

#include "global_quantities_kerr.h"
#include "global_quantities_main.h"
#include "global_quantities_osc.h"
#include "physical_quantities.h"
#include "macros_2.h"
//#include "macros_1.h"


int analytic_r_costh(double*,double *, double *);
int analytic_derivatives(double *,double *,double *,double *,double *,double *,double *,double *,double *, double *, double *, double *, double *, double *,double *, double *, double *, double *,double *,double *,double *,double *,double , double );
int compute_frequencies(double *, double *, double *, double *, double *,double *, double *, double *, double *);
int dchi_dpsi_dphi_dt(double *, double, double,double *,double *,double *,double *);
int pez_to_ELCQ(char, double, double, double,double, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int PN_ELzQdot(double *, double *, double *, double *, double, double,double,double,double,double);
int turning_points(double, double, double, double, double, double, double*,  double *, double *, double *, double *);
int turning_points_derivatives(double*, double* ,double* ,double* ,double* , double *, double* ,double*, double*,double*);

double *Realvector(const long nl, const long nh);
void free_Realvector(double *v, const long nl, const long nh);


int AA_RHS(double *x,double *dxdt)
{
  double *drpd, *drad,*dr3d,*dr4d; // Vector of radial Turning points derivatives
  double *dzmd,*dzpd;              // Vector of polar Turning points derivatives
  double *decd, *dpd, *dthind ;   // Orbital elements evolution equation
  double *dOmRd, *dOmTHd;          // Orbital (Mino) frequencies derivatives
  double drpdE, dradE, drpdLz, dradLz, drpdQ, dradQ; // "
  double dr3dE, dr4dE, dr3dLz, dr4dLz, dr3dQ, dr4dQ; // "
  double dzmdE, dzpdE, dzmdLz, dzpdLz, dzmdQ, dzpdQ; // "
  double drdra, drdrp, drdr3, drdr4,dcosthdzp,dcosthdzm; // Radial and polar derivatives wrt  turning points
  double omE2,beta,sigma, delta;            // Auxiliar variables
  double dVrdr, dVthdcth;                   // Derivatives of the radial and polar potentials
  double drdw, dcosthdw;                    // Radial and plolar derivatives wrt their fundamental frequencies
  double d2rdlE,d2rdlLz,d2rdlQ;             // Derivatives of the radial "geodesic eq." wrt constants of the motion
  double d2costhdlE,d2costhdlQ,d2costhdlLz; // Derivatives of the polar "geodesic eq" wrt constants of the motion
  double r, costh;
  double dpdE, dpdLz, dpdQ,decdE, decdLz, decdQ, dthindE, dthindLz, dthindQ;
  double dpsidt, dchidt;
  int Ia;
  
  Ia = 3;

//--- Allocating memory for vector derivatives
  drpd = Realvector(1, Ia);
  drad = Realvector(1, Ia);
  dr3d = Realvector(1, Ia);
  dr4d = Realvector(1, Ia);
  
  dzmd = Realvector(1, Ia);
  dzpd = Realvector(1, Ia);
    
  decd = Realvector(1, Ia);
  dpd  = Realvector(1, Ia);
  dthind = Realvector(1, Ia);
   
  dOmRd = Realvector(1, Ia);
  dOmTHd = Realvector(1, Ia);
    
/*-----Computing new Constants of Motion and "main" turning points-----*/
 pez_to_ELCQ(FLAG_CircOrb, x[3], x[4], x[5],SPIN,&INCLINATION, &R_peri, &R_apo, &THETA_min, &Z_minus, &ENERGY, &LZ, &C_CONSTANT,&Q_CONSTANT);
 
  omE2 = 1.0-SQR(ENERGY);
  beta = SQR(SPIN)*omE2;

/*-----Computing the rest of the turning points -----*/
  turning_points(x[3], x[4], ENERGY, LZ,Q_CONSTANT,Z_minus, &R_3, &R_4, &P3_p, &P4_p, &Z_plus);
   
/*----Computing new fundamental frequencies (associated with Mino time) and initial arguments of the Elliptic functions ---*/
  compute_frequencies(&OMEGA_r, &OMEGA_th, &OMEGA_t, &Kcom_r, &Kcom_th, &Ecom_r, &Ecom_th, &KR, &KTH);
       
/*---New Radial and polar coordinates parameterized in function of the action angle variables---*/
  analytic_r_costh(x, &r, &costh);
  
  sigma = r*r + SPIN*SPIN*costh*costh;
  delta = r*r + SPIN*SPIN*costh*costh - 2.0*r;
    
  
/*--- Computing the derivatives of the orbital parameters and the Constats of motion----*/
  turning_points_derivatives(drpd,drad,dr3d,dr4d,dzmd,dzpd,dpd,decd,dthind,x);

/*----- Computing first and second derivatives of the radial and polar coordinates -----*/
  analytic_derivatives(&d2rdlE, &d2rdlLz, &d2rdlQ, &d2costhdlE, &d2costhdlQ, &d2costhdlLz,
                       drpd, drad, dr3d, dr4d, dzmd, dzpd,dOmRd,dOmTHd,
                       &drdra, &drdrp, &drdr3, &drdr4, &dcosthdzp, &dcosthdzm,
                       &drdw, &dcosthdw,
                       x[6], x[7]);
  
  
  drpdE = drpd[0]; dradE = drad[0]; drpdLz = drpd[1]; dradLz = drad[1]; drpdQ = drpd[2]; dradQ = drad[2];
  dr3dE = dr3d[0]; dr4dE = dr4d[0]; dr3dLz = dr3d[1]; dr4dLz = dr4d[1]; dr3dQ = dr3d[2]; dr4dQ = dr4d[2];
  dzmdE = dzmd[0]; dzpdE = dzpd[0]; dzmdLz = dzmd[1]; dzpdLz = dzpd[1]; dzmdQ = dzmd[2]; dzpdQ = dzpd[2];
  
  decdE = decd[0]; decdLz = decd[1]; decdQ = decd[2];
  dpdE = dpd[0];   dpdLz = dpd[1];   dpdQ = dpd[2];
  dthindE = dthind[0]; dthindLz = dthind[1]; dthindQ = dthind[2]; // Z_minus check

/*----- Derivatives of the radial and polar "potentials" -----*/
  dVrdr = omE2*( ((R_apo - r)-(r - R_peri))*(r - R_3)*(r-R_4)
                +(R_apo - r)*(r - R_peri)*( (r-R_4) + (r - R_3 )));
  
  dVthdcth = -beta*2.*costh*( (Z_plus - SQR(costh)) + (Z_minus - SQR(costh)) );
  
  
/*----- Computing velocities d(Chi,Psi,Phi, dTau)/dt-----*/
  dchi_dpsi_dphi_dt(x, r, costh, &dpsidt, &dchidt,&dxdt[1], &dxdt[2]);
  

  dxdt[8] = dxdt[2]/sigma; // dlambdadt  CHECK TIME ODEs
    
/*---- Evolution of the constants of motion using PN self-force [instant changes in MBH units] ---*/
  PN_ELzQdot(&ENERGY_dot, &LZ_dot, &Q_CONSTANT_dot,x,r, costh, dpsidt, dchidt,dxdt[1], dxdt[2]);
  

/*--- Evolution equations for the orbital elements---*/
  //---dp_dt
  dxdt[3] = (dpdE*ENERGY_dot + dpdLz*LZ_dot + dpdQ*Q_CONSTANT_dot);
  //---dec_dt
  dxdt[4] = (decdE*ENERGY_dot+ decdLz*LZ_dot+ decdQ*Q_CONSTANT_dot);
  //---dthinc_dt
  dxdt[5] = (dthindE*ENERGY_dot+dthindLz*LZ_dot+dthindQ*Q_CONSTANT_dot);
  
  //---dwr_dt
  dxdt[6] = (dOmRd[0]*ENERGY_dot + dOmRd[1]*LZ_dot + dOmRd[2]*Q_CONSTANT_dot)*x[8] + OMEGA_r*dxdt[8];
  
  //---dwth_dt
  dxdt[7] = (dOmTHd[0]*ENERGY_dot + dOmTHd[1]*LZ_dot + dOmTHd[2]*Q_CONSTANT_dot)*x[8] + OMEGA_th*dxdt[8];
  
 
  // printf("ENERGY_dot=%4.6e, LZ_dot=%4.6e, Q_CONSTANT_dot=%4.6e\n",ENERGY_dot, LZ_dot, Q_CONSTANT_dot);
  
  // printf("dxdt[6]=%4.6e dxdt[7] =%4.6e\n",dxdt[6],dxdt[7]);
  

  
/*
  // ==== OSCULATING EVOLUTION EQUATIONS ===//
  
  double TOL;
  TOL =1.0E-6;
  
   
  if (fabs(drdw) > TOL)
  {
    dxdt[6] =( drdr3*(dr3dE*ENERGY_dot + dr3dQ*Q_CONSTANT_dot) + drdr4*(dr4dE*ENERGY_dot + dr4dQ*Q_CONSTANT_dot));
    dxdt[6] /= -drdw;
  }
  else
  {
    dxdt[6]  = OMEGA_r*( -PNFcomp[1] + 2.0*(d2rdlE*ENERGY_dot + d2rdlQ*Q_CONSTANT_dot)*dxdt[8]);
    dxdt[6]/= dVrdr;
  }
  
  if (fabs(dcosthdw) > TOL)
  {
     dxdt[7]  = (dcosthdzp*(dzpdE*ENERGY_dot + dzpdLz*LZ_dot + dzpdQ*Q_CONSTANT_dot));
     dxdt[7] /= -dcosthdw;
  }
  else
  {
     dxdt[7]  = OMEGA_th*( -PNFcomp[2] + 2.0*(d2costhdlE*ENERGY_dot + d2costhdlLz*LZ_dot + d2costhdlQ*Q_CONSTANT_dot)*dxdt[8] );
     dxdt[7] /= dVthdcth;
  }

*/

//--- Asigning global variables
  dRdt = drdw*OMEGA_r*dxdt[8];
  dCOSTHdt = dcosthdw*OMEGA_th*dxdt[8];

//--- Freeing memory for vector derivatives
  free_Realvector(drad, 1, Ia);
  free_Realvector(drpd, 1, Ia);
  free_Realvector(dr3d, 1, Ia);
  free_Realvector(dr4d, 1, Ia);
  free_Realvector(dzmd, 1, Ia);
  free_Realvector(dzpd, 1, Ia);
  free_Realvector(decd, 1, Ia);
  free_Realvector(dpd, 1, Ia);
  free_Realvector(dthind, 1, Ia);
  free_Realvector(dOmRd,1,Ia);
  free_Realvector(dOmTHd,1,Ia);
/*----- This is the end -----*/
  return 0;
}


