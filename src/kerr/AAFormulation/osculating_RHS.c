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
 Last Update: 29/07/2013                                      
 =============================================================================*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"

#include "physical_quantities.h"

//#include "macros_1.h"
#include "macros_2.h" 

int analytic_r_costh(double*,double *, double *);
int analytic_derivatives(double*, double*, double*, double*, double*, double*,double*,double*, double*, double*, double*, double*,double*,double*,double, double);
int compute_frequencies(double *, double *, double *, double *, double *,double *, double *, double *, double *);
int pez_to_ELCQ(char, double, double, double,double, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int PN_ELzQdot(double *, double *, double *, double *, double *, double, double,double,double,double,double);
int turning_points(double, double, double, double, double, double, double*, double *,  double *, double *, double *);
int turning_points_derivatives(double*, double* ,double* ,double* ,double* ,double *, double *, double* ,double*, double* );

int dchi_dpsi_dphi_dt(double *, double *, double, double,double *,double *,double *,double *);

double *Realvector(const long nl, const long nh);
void free_Realvector(double *v, const long nl, const long nh);


int osculating_RHS(double *x,double *dxdt)
{
  double *drpd, *drad,*dr3d,*dr4d; // Vector of radial Turning points derivatives
  double *dzmd,*dzpd;              // Vector of polar Turning points derivatives
  double *decd, *dpd, *diotad;      // Orbital elements evolution equation
  
  double dWrodt, dWTHodt, Wr, Wth; // Angle variables and derivatives
  double drpdE, dradE, drpdLz, dradLz, drpdQ, dradQ; // "
  double dr3dE, dr4dE, dr3dLz, dr4dLz, dr3dQ, dr4dQ; // "
  double dzmdE, dzpdE, dzmdLz, dzpdLz, dzmdQ, dzpdQ; // "
  double drdra, drdrp, drdr3, drdr4,dcosthdzp,dcosthdzm; // Radial and polar derivatives wrt  turning points
  double omE2,beta, SomE2, sigma;                  // Auxiliar variables
  double dVrdr, dVthdcth;                   // Derivatives of the radial and polar potentials
  double drdw, dcosthdw;                    // Radial and plolar derivatives wrt their fundamental frequencies
  double d2rdlE,d2rdlLz,d2rdlQ;             // Derivatives of the radial "geodesic eq." wrt constants of the motion
  double d2costhdlE,d2costhdlQ,d2costhdlLz; // Derivatives of the polar "geodesic eq" wrt constants of the motion
  double TOL;                               // Tolerance set to choose the oscullating condition to be applied.
  double r, costh;
  double dpdE, dpdLz, dpdQ,decdE, decdLz, decdQ, diotadE, diotadLz, diotadQ;
  double dpsidt, dchidt,dphidt, dtaudt;
  int i,Ia;
  
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
  diotad = Realvector(1, Ia);
   

/*-----Computing new Constants of Motion and "main" turning points-----*/
  pez_to_ELCQ(FLAG_circular_orbit, x[3], x[4], THETA_inc,SPIN,
              &x[5], &R_peri, &R_apo, &THETA_min, &Z_minus,
              &ENERGY, &LZ, &C_CONSTANT, &Q_CONSTANT);
    
  omE2 = 1.0-SQR(ENERGY);
  beta = SQR(SPIN)*omE2;
  SomE2 = SQR(omE2);
    
/*-----Computing the rest of the turning points -----*/
  turning_points(x[3], x[4], ENERGY, LZ,Q_CONSTANT,Z_minus, &R_3, &R_4, &P3_p, &P4_p, &Z_plus);

/*--- Computing the derivatives of the orbital parameters and the Constats of motion----*/
  turning_points_derivatives(drpd,drad,dr3d,dr4d,dzmd,dzpd,dpd,decd,diotad,x);
  
/*----Computing new fundamental frequencies (associated with Mino time) and initial arguments of the Elliptic functions ---*/
  compute_frequencies(&OMEGA_r, &OMEGA_th, &OMEGA_t, &Kcom_r, &Kcom_th, &Ecom_r, &Ecom_th, &KR, &KTH);

/*---Radial and polar coordinates parameterized in function of the action angle variables---*/
  analytic_r_costh(x, &r, &costh);
  
/*----- Computing velocities d(Chi,Psi,Phi, dTau)/dt-----*/
  dchi_dpsi_dphi_dt(dxdt,x, r, costh, &dpsidt, &dchidt,&dphidt, &dtaudt);
  
/*---- Evolution of the constants of motion using PN self-force [instant changes in MBH units] ---*/
  PN_ELzQdot(&ENERGY_dot, &LZ_dot, &Q_CONSTANT_dot,dxdt,x,r, costh, dpsidt, dchidt,dphidt, dtaudt);
   
/*----- Derivatives of the radial and polar "potentials" -----*/
  dVrdr = omE2*( ((R_apo - r)-(r - R_peri))*(r - R_3)*(r-R_4)
                 +(R_apo - r)*(r - R_peri)*( (r-R_4) + (r - R_3 )));
    
  dVthdcth = beta*2.*costh*(1.- SQR(costh))*( (Z_plus - SQR(costh)) + (Z_minus - SQR(costh)) );
  
/*----- Computing first and second derivatives of the radial and polar coordinates -----*/
  analytic_derivatives(&d2rdlE, &d2rdlLz, &d2rdlQ, &d2costhdlE, &d2costhdlQ, &d2costhdlLz,
                         drpd, drad, dr3d, dr4d, dzmd, dzpd,
                         &drdw, &dcosthdw,
                         x[6], x[7]);

 
/*===========================
   Osculating Conditions
============================*/
  // (1) First set of oscullating conditions: Singular at ∂r/∂w = 0
  // NOTE: On resonace the turnining points for r and theta don't coincide!!
  // (2) Second set of oscullating conditions: Singular at ∂V/∂q_µ = 0 NOTE: CHECK THE 2 FACTOR
  
  drpdE = drpd[0]; dradE = drad[0]; drpdLz = drpd[1]; dradLz = drad[1]; drpdQ = drpd[2]; dradQ = drad[2];
  dr3dE = dr3d[0]; dr4dE = dr4d[0]; dr3dLz = dr3d[1]; dr4dLz = dr4d[1]; dr3dQ = dr3d[2]; dr4dQ = dr4d[2];
  dzmdE = dzmd[0]; dzpdE = dzpd[0]; dzmdLz = dzmd[1]; dzpdLz = dzpd[1]; dzmdQ = dzmd[2]; dzpdQ = dzpd[2];

  decdE = decd[0]; decdLz = decd[1]; decdQ = decd[2];
  dpdE = dpd[0];   dpdLz = dpd[1];   dpdQ = dpd[2];
  diotadE = diotad[0]; diotadLz = diotad[1]; diotadQ = diotad[2]; // Z_minus check
    
  printf("Edt=%4.6e, Lz_dt=%4.6e, Q_dt=%4.6e \n", ENERGY_dot, LZ_dot, Q_CONSTANT_dot);
  printf(" pdot=%4.6e, edot=%4.6e, iotadot=%4.6e \n\n", dxdt[3], dxdt[4], dxdt[5]);
  
//--- Evolution equations for the orbital elements
    //---p_dot
    dxdt[3] = (dpdE*ENERGY_dot + dpdLz*LZ_dot + dpdQ*Q_CONSTANT_dot);
    //---ec_dot
    dxdt[4] = (decdE*ENERGY_dot+ decdLz*LZ_dot+ decdQ*Q_CONSTANT_dot);
    //---iota_dot
    dxdt[5] = -(diotadE*ENERGY_dot+diotadLz*LZ_dot+diotadQ*Q_CONSTANT_dot)*sin(INCLINATION);
  
    printf("dpdE=%4.6e, dpdLz=%4.6e, dpdQ=%4.6e \n",dpdE, dpdLz, dpdQ);
    printf("decdE=%4.6e, drcdLz=%4.6e, decdQ=%4.6e\n",decdE, decdLz, decdQ);
    printf("didE=%4.6e, didLz=%4.6e, didQ=%4.6e \n",diotadE, diotadLz, diotadQ);
  
  
  
  printf(" drdw=%4.6e, dVrdr=%4.6e dcosthdwdr=%4.6e Vthdcthr=%4.6e\n ",drdw, dVrdr,dcosthdw, dVthdcth);
  printf("PNFcomp[1]=%4.6e, PNFcomp[2]=%4.6e,PNFcomp[3]=%4.6e\n\n",PNFcomp[1], PNFcomp[2],PNFcomp[3]);
  
  TOL = 1.0E-9;
  
  if (fabs(drdw) > TOL)
  {
    dxdt[6] = ((drdra*dradE  + drdrp*drpdE  + drdr3*dr3dE  + drdr4*dr4dE)*ENERGY_dot
               +(drdra*dradLz + drdrp*drpdLz + drdr3*dr3dLz + drdr4*dr4dLz)*LZ_dot
               +(drdra*dradQ  + drdrp*drpdQ  + drdr3*dr3dQ  + drdr4*dr4dQ)*Q_CONSTANT_dot
             );
    dxdt[6] /= drdw; //check this
    
     printf("wr 1\n");
  }
  else
  {
    dxdt[6] = OMEGA_r*( -PNFcomp[1] + 2.0*(d2rdlE*ENERGY_dot + d2rdlLz*LZ_dot + d2rdlQ*Q_CONSTANT_dot));
    dxdt[6] /= dVrdr;
    

     printf("wr 2\n");
  }
  
    
  if (fabs(dcosthdw) > TOL)
  {
    dxdt[7] = ( dcosthdzm*(dzmdE*ENERGY_dot + dzmdLz*LZ_dot + dzmdQ*Q_CONSTANT_dot)
               +dcosthdzp*(dzpdE*ENERGY_dot + dzpdLz*LZ_dot + dzpdQ*Q_CONSTANT_dot));
    dxdt[7] /= dcosthdw;
    
     printf("wth 1\n");
  }
  else
  {
    dxdt[7] = OMEGA_th*( -PNFcomp[2] + 2.0*(d2costhdlE*ENERGY_dot + d2costhdlLz*LZ_dot + d2costhdlQ*Q_CONSTANT_dot)*dxdt[4]/sigma);
    dxdt[7] /= dVthdcth;

    printf("wth 2\n");
  }
  
   
   //  printf(" \n ");

    
    
    
  
  //printf("P_p=%4.6e,ECCENTRICITY=%4.6e,THETA_inc=%4.6e, R_p_o=%4.6e, COSTH_p_o=%4.6e, ENERGY=%4.6e, LZ=%4.6e, Q_CONSTANT=%4.6e,INCLINATION=%4.6e, PSI_NEW=%4.6e, CHI_NEW=%4.6e\n\n", P_p,ECCENTRICITY,THETA_inc, R_p_o, COSTH_p_o, ENERGY, LZ, Q_CONSTANT,INCLINATION, PSI_NEW, CHI_NEW);
  
  
 //  printf("ENERGY_dot=%4.6e, LZ_dot=%4.6e,Q_CONSTANT_dot=%4.6e\n\n",ENERGY_dot, LZ_dot,Q_CONSTANT_dot, );

    
//--- Freeing memory for vector derivatives
  free_Realvector(drad, 1, Ia);
  free_Realvector(drpd, 1, Ia);
  free_Realvector(dr3d, 1, Ia);
  free_Realvector(dr4d, 1, Ia);
  free_Realvector(dzmd, 1, Ia);
  free_Realvector(dzpd, 1, Ia);
  free_Realvector(decd, 1, Ia);
  free_Realvector(dpd, 1, Ia);
  free_Realvector(diotad, 1, Ia);
  
/*----- This is the end -----*/
  return 0;
}


