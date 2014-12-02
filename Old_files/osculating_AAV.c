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
 find angle variables w_r = omega_r*t=Omega_r*λ, w_theta = omega_theta*t
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
                                            Created by Priscilla on 12/05/2013
 -------------------------------------------------------------------------------
 Last Update: 08/04/2013                                        Modified by PCM
 =============================================================================*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "global_quantities_resonances.h"

#include "physical_quantities.h"

//#include "macros_1.h"
#include "macros_2.h"

void JacobianEllipticFunctions(double, double, double *, double * , double * );

int osculating_RHS(double *dwrdl, double *dwthdl, double *dphidl, double w_r, double w_th)
{
    double dr3dE, dr4dE, dr3dQ, dr4dQ;                // Derivatives of the r3 and r4 turning poits
    double dkrdr3, dkrdr4, dkthdzp;                   // Derivatives of the argument of the eliptical function
    double dzpdE, dzpdQ,dzpdL;                        // Derivatives of the polar turning pints
    double one_minus_E2, Sone_minus_E2;               // Auxiliar varialbles
    double dElKdkr,dElKdkth;                          // Derivatives elliptical functions
    double sn_r, sn_th;                               // Jacobi's elliptic functions
    double cn_r, cn_th;                               // Jacobi's elliptic functions
    double dn_r, dn_th;                               // Jacobi's elliptic functions
    double varphi_r,varphi_th, dvp_dzp;               // Argument Jacobi's elliptic functions
    double dvp_r3, dvp_r4;                            // Derivatives argument Jacobi elliptical function
    double A, B, dAdr3, dAdr4,dBdr3, dBdr4;           // denominator and nominator of radial coordinates and derivatives
    double drdr3,drdr4;                               // Radial derivative wrt the r3 and r4 turning poins
    double dcosthdzp;                                 // Polar derivative wrt the z+ turning point
    double drdw, dcosthdw;                            // Radial and plolar derivatives wrt their fundamental frequencies
    double denom, r3_root;                            // Aux variables
    double r2pa2,delta,SQr2pa2,gamma;                 //   "
    double denom_phi,numer_phi;                       //   "
    double drdl, dcosthdl;                            // Geodesic equations for the BL coordinates r and costh in function of w_r, w_th respectively
    double rp, cosTh,theta;                           // Radial and polar BL coordinates
    double g1,N,D, dNdr3, dNdr4, dDdE, dDdQ;          // Aux variables
    double dOmegaRdE,dOmegaRdQ;                       // Derivatives of the fundamental frequencies wrt constants of the motion
    double dOmegaTHdE, dOmegaTHdL, dOmegaTHdQ;        //  "
    double epsilon_o;                                 // Aux variable
    double d2rdlE,d2rdlQ;                             // Derivatives of the radial "geodesic eq." wrt constants of the motion
    double d2costhdlE,d2costhdlQ,d2costhdlL;          // Derivatives of the polar "geodesic eq" wrt constants of the motion
    double dVrdr, dVthdcth;                           // Derivatives of the radial and polar potentials
    double TOL;                                       // Tolerance set to choose the oscullating condition to be applied.
    int sign_r, sign_th;                              // Aux variables

/*----- Computind dPhi_dl -----*/
    rp    = GV_R_p_o;
    cosTh = GV_CosTh_p_o;
    
    theta = acos(cosTh);
    
    r2pa2= SQR(rp) + GV_Spin2;                        // r^2 + a^2
    SQr2pa2 = SQR(r2pa2);                             // (r^2+a^2)^2
    delta =  r2pa2 -2*rp;
    gamma = GV_E*(SQr2pa2/delta-GV_Spin2) - 2.0*GV_Spin*rp*GV_Lz/delta; // gamma = E [ ( r^2 + a^2 )^2 / Delta - a^2 ] - 2 M a r Lz / Delta
    
    numer_phi = ( (1.0/(1.0-SQR(cosTh)) - GV_Spin2/delta)*GV_Lz + 2.0*GV_Spin*GV_E*rp/delta );
    denom_phi = sqrt((GV_Spin2*(1.0 - SQR(GV_E)))*(GV_Z_plus-SQR(cosTh)));
   
    //--Azimutal geodesic equation
    *dphidl = numer_phi/denom_phi;
    
/*--- Osculating evolution for the radial and polar variables---*/    
    
//---Computing the derivatives of the Elliptic functions' argument 
    dkrdr3 =  (GV_R_apo - GV_R_peri)*(GV_R_apo-GV_R_4);
    dkrdr3 /= 2.0*GV_kr*(GV_R_peri-GV_R_4)*SQR(GV_R_apo - GV_R_3);
    
    dkrdr4 = (GV_R_apo - GV_R_peri)*(GV_R_3-GV_R_peri);
    dkrdr4 /= 2.0*GV_kr*(GV_R_apo - GV_R_3)*SQR(GV_R_peri-GV_R_4);
    
    dkthdzp =-GV_Z_minus/(2.0*GV_kth*SQR(GV_Z_plus));
    
//---Computing the partial derivatives of the turning points wrt constants of motion
    one_minus_E2 = (1.0-SQR(GV_E));
    Sone_minus_E2 = SQR(one_minus_E2);
    epsilon_o = GV_Spin2*(1.0-SQR(GV_E))/SQR(GV_Lz);
    
    r3_root = sqrt(SQR(2.0/one_minus_E2 - (GV_R_peri+GV_R_apo)) -4.0*SQR(GV_Spin)*GV_Q/(one_minus_E2*GV_R_apo*GV_R_peri));
    
    dr3dE = 2.0*GV_E*(1.0 + (1.0/one_minus_E2 - SQR(GV_Spin)*GV_Q/(GV_R_peri*GV_R_apo))/r3_root  );
    dr3dE /= Sone_minus_E2;
    
    dr4dE = SQR(GV_Spin)*GV_Q*(2.0*GV_E*GV_R_3/Sone_minus_E2-dr3dE/one_minus_E2);
    dr4dE /= GV_R_peri*GV_R_apo*SQR(GV_R_3);
    
    dr3dQ = -SQR(GV_Spin);
    dr3dQ /= one_minus_E2*GV_R_peri*GV_R_apo*r3_root;
    
    dr4dQ = SQR(GV_Spin)*(GV_R_3 -GV_Q*dr3dQ);
    dr4dQ /= one_minus_E2*GV_R_peri*GV_R_apo*SQR(GV_R_3);
    
    dzpdE = 2.0*SQR(GV_Spin2)*GV_Q*GV_E;
    dzpdE /=GV_Z_minus*SQR(epsilon_o*GV_Lz)*SQR(GV_Lz);
    
    dzpdQ = 1.0/(SQR(GV_Spin)*GV_Z_minus*epsilon_o);
    
    dzpdL = -2.0*GV_Q*(epsilon_o*GV_Lz -GV_Spin2*one_minus_E2/GV_Lz);
    dzpdL = GV_Z_minus*SQR(epsilon_o*GV_Lz)*SQR(GV_Lz);
    
    // Sanity checks
    /*double r3= 0.5*(2/one_minus_E2-(GV_R_peri+GV_R_apo) +r3_root);
     double r4  = SQR(GV_Spin)*GV_Q/(one_minus_E2*GV_R_apo*GV_R_peri*r3);
     double r3_p = GV_P_3/(1-GV_Ecc);
     double r4_p = GV_P_4/(1+GV_Ecc);
     printf("r3=%4.6e, r3_p=%4.6e, R3=%4.6e, r4= %4.6e, r4_p= %4.6e, R4= %4.6e\n", r3, r3_p, GV_R_3, r4, r4_p, GV_R_4);
     */
    
    
/*---- Derivatives of the Complete elliptical integrals K(k_r), K(k_th)  ---*/
    dElKdkr = GV_Ecomp_r/(GV_kr*(  1.0-SQR(GV_kr))) - GV_Kcomp_r/GV_kr;
    dElKdkth = GV_Ecomp_th/(GV_kth*(1.0-SQR(GV_kth)))- GV_Kcomp_th/GV_kth;
    
/*--- Arguments of the Jacobi elliptic functions sn(varphi,k), cn(varphi,k), dn(varphi,k) ----*/
    sign_r =0;
    sign_th=0;
    
    // Radial argument
    if( (w_r>= 0)&&(w_r<=PI) )
    {
        varphi_r = w_r*GV_Kcomp_r/PI;
        sign_r = 1;
    }
    else if( (w_r>= PI)&&(w_r<=2*PI) )
    {
        varphi_r = (2*PI-w_r)*GV_Kcomp_r/PI;
        sign_r = -1;
    }
    else
        printf("w_r out of bounds!!\n");
    
    // Polar argument
    if( (w_th>= 0)&&(w_th<=PI/2) )
    {
        varphi_th = w_th*2*GV_Kcomp_th/PI;
        sign_th = 1;
    }
    else if( (w_th>= PI/2)&&(w_th<=3*PI/2) )
    {
        varphi_th = (PI-w_th)*2*GV_Kcomp_th/PI;
        sign_th =-1;
    }
    else if( (w_th>= 3*PI/2)&&(w_th<=PI) )
    {
        varphi_th = (w_th-2*PI)*2*GV_Kcomp_th/PI;
        sign_th = 1;
    }
    else
        printf("w_th out of bounds!!\n");
    
/*---Jacobi elliptic functions----*/
    JacobianEllipticFunctions(varphi_r, SQR(GV_kr),  &sn_r,  &cn_r ,  &dn_r);
    JacobianEllipticFunctions(varphi_th,SQR(GV_kth), &sn_th, &cn_th , &dn_th);
    
/*---- Derivatives of the arguments of the Jacobi elliptical function sn(varphi,k), cn(varphi,k), dn(varphi,k) ---*/
    if( (w_r>= 0)&&(w_r<=PI) )
    {
        dvp_r3 = w_r*dElKdkr*dkrdr3;
        dvp_r3 /= PI;
        
        dvp_r4 = w_r*dElKdkr*dkrdr4;
        dvp_r4 /= PI;
    }
    else if( (w_r>= PI)&&(w_r<=2*PI) )
    {
        dvp_r3 = (2.0*PI-w_r)*dElKdkr*dkrdr3;
        dvp_r3 /= PI;
        
        dvp_r4 = (2.0*PI-w_r)*dElKdkr*dkrdr4;
        dvp_r4 /= PI;
    }
    
    if( (w_th>= 0)&&(w_th<=PI/2) )
        dvp_dzp = w_th*2.0*dElKdkth*dkthdzp/PI;
    else if( (w_th>= PI/2)&&(w_th<=3*PI/2) )
        dvp_dzp = (PI-w_th)*2.0*dElKdkth*dkthdzp/PI;
    else if( (w_th>= 3*PI/2)&&(w_th<=PI) )
        dvp_dzp = (w_th-2*PI)*2*dElKdkth*dkthdzp/PI;
    
/*---- Partial derivatives of r(w_b)  [r = A/B ] wrt the turning points r3, r4----*/
    A = GV_R_3*(GV_R_apo - GV_R_peri)*SQR(sn_r) - GV_R_peri*(GV_R_apo - GV_R_3);
    B = (GV_R_apo - GV_R_peri)*SQR(sn_r) -(GV_R_apo - GV_R_3);
    
    dAdr3 = (GV_R_apo - GV_R_peri)*(SQR(sn_r) + GV_R_3*2.0*sn_r*cn_r*dn_r*dvp_r3) + GV_R_peri;
    dAdr4 = GV_R_3*(GV_R_apo - GV_R_peri)*2.0*sn_r*cn_r*dn_r*dvp_r4;
    
    dBdr3 = (GV_R_apo - GV_R_peri)*2.0*sn_r*cn_r*dn_r*dvp_r3 + 1.0;
    dBdr4 = (GV_R_apo - GV_R_peri)*2.0*sn_r*cn_r*dn_r*dvp_r4;
    
    drdr3 = B*dAdr3-A*dBdr3;
    drdr3 /=SQR(B);
    
    drdr4 = B*dAdr4-A*dBdr4;
    drdr4 /=SQR(B);
    
/*---- Partial derivatives of cos(theta)(w_th) wrt turning points----*/
    dcosthdzp = sqrt(GV_Z_minus)*cn_th*dn_th*dvp_dzp;
    
/*---- Computing the second derivatives of the geodesic equations [dr/dλ = N*D/B^2, dcos(theta)[w_th]/dλ] wrt the constant of motion---*/
    
    g1 = (GV_R_peri-GV_R_3)*(GV_R_apo-GV_R_3)*(GV_R_apo-GV_R_peri);
    N  = sn_r*cn_r*dn_r*g1;
    D  = GV_Kcomp_r*GV_Omega_r/PI;
    
    dNdr3 = (SQR(cn_r*dn_r) - SQR(sn_r)*( SQR(dn_r) + SQR(GV_kr*cn_r) ) )*g1*dvp_r3 - sn_r*cn_r*dn_r*(GV_R_apo+GV_R_peri)*(GV_R_apo-GV_R_peri);
    dNdr4 = (SQR(cn_r*dn_r) - SQR(sn_r)*( SQR(dn_r) + SQR(GV_kr*cn_r) ) )*g1*dvp_r4;
    
    dOmegaRdE = -SQR(PI)*(2.0*GV_E*(GV_R_apo-GV_R_3)*(GV_R_peri-GV_R_4) +(1.0-SQR(GV_E))*( (GV_R_peri-GV_R_4)*dr3dE + (GV_R_apo-GV_R_3)*dr4dE));
    dOmegaRdE /= GV_Omega_r*8.0*SQR(GV_Kcomp_r);
    dOmegaRdE += -GV_Omega_r*dElKdkr*(dkrdr3*dr3dE + dkrdr4*dr4dE)/GV_Kcomp_r;
    
    dOmegaRdQ = -SQR(PI)*( (1.-SQR(GV_E))*( (GV_R_peri-GV_R_4)*dr3dQ + (GV_R_apo-GV_R_3)*dr4dQ));
    dOmegaRdQ /= GV_Omega_r*8.0*SQR(GV_Kcomp_r);
    dOmegaRdQ += -GV_Omega_r*dElKdkr*(dkrdr3*dr3dQ + dkrdr4*dr4dQ)/GV_Kcomp_r;
    
    dOmegaTHdE =  PI*GV_Lz*(-GV_Z_plus*(2.0*GV_Spin2*GV_E)/SQR(GV_Lz) + epsilon_o*dzpdE );
    dOmegaTHdE /= 4.0*GV_Kcomp_th*sqrt(epsilon_o*GV_Z_plus);
    dOmegaTHdE += -GV_Omega_th*dElKdkth*dkthdzp*dzpdE/GV_Kcomp_th;
    
    dOmegaTHdQ =  PI*GV_Lz*epsilon_o*dzpdQ;
    dOmegaTHdQ /= 4.0*GV_Kcomp_th*sqrt(epsilon_o*GV_Z_plus);
    dOmegaTHdQ += -GV_Omega_th*dElKdkth*dkthdzp*dzpdQ/GV_Kcomp_th;
    
    dOmegaTHdL =  PI*GV_Lz*(-GV_Z_plus*(2.0*GV_Spin2*GV_E)/SQR(GV_Lz) + epsilon_o*dzpdL );
    dOmegaTHdL /= 4.0*GV_Kcomp_th*sqrt(epsilon_o*GV_Z_plus);
    dOmegaTHdL += (-GV_Omega_th*dElKdkth*dkthdzp*dzpdL +0.5*PI*sqrt(epsilon_o*GV_Z_plus))/GV_Kcomp_th;
    
    dDdE = dElKdkr*(dkrdr3*dr3dE + dkrdr4*dr4dE)*GV_Omega_r/PI + GV_Ecomp_r*dOmegaRdE;
    dDdQ = dElKdkr*(dkrdr3*dr3dQ + dkrdr4*dr4dQ)*GV_Omega_r/PI + GV_Ecomp_r*dOmegaRdQ;
    
    
    d2rdlE  = (2.0*sign_r*( (dNdr3*dr3dE + dNdr4*dr4dE)*SQR(B) - 2.0*N*B*(dBdr3*dr3dE + dBdr4*dr4dE))*D)/SQR(B*B);
    d2rdlE += N*(dElKdkr*(dkrdr3*dr3dE + dkrdr4*dr4dE )*GV_Omega_r + GV_Kcomp_r*dOmegaRdE)/(PI*SQR(B));
    
    d2rdlQ  = (2.0*sign_r*( (dNdr3*dr3dQ + dNdr4*dr4dQ)*SQR(B) - 2.0*N*B*(dBdr3*dr3dQ + dBdr4*dr4dQ))*D)/SQR(B*B);
    d2rdlQ += N*(dElKdkr*(dkrdr3*dr3dQ + dkrdr4*dr4dQ )*GV_Omega_r + GV_Kcomp_r*dOmegaRdQ)/(PI*SQR(B));

    
    d2costhdlE = 2.0*sign_th*sqrt(GV_Z_minus)*( (-(sn_th*SQR(dn_r) + SQR(GV_kth*cn_th)*sn_th)*GV_Kcomp_th*GV_Omega_th*dvp_dzp*dzpdE
                                              +cn_th*dn_th*(dElKdkth*dkthdzp*dzpdE*GV_Omega_th + GV_Kcomp_th*dOmegaTHdE)))/PI;
    d2costhdlL = 2.0*sign_th*sqrt(GV_Z_minus)*(-(sn_th*SQR(dn_r) + SQR(GV_kth*cn_th)*sn_th)*GV_Kcomp_th*GV_Omega_th*dvp_dzp*dzpdL
                                              +cn_th*dn_th*(dElKdkth*dkthdzp*dzpdL*GV_Omega_th + GV_Kcomp_th*dOmegaTHdL))/PI;
    d2costhdlQ = 2.0*sign_th*sqrt(GV_Z_minus)*(-(sn_th*SQR(dn_r) + SQR(GV_kth*cn_th)*sn_th)*GV_Kcomp_th*GV_Omega_th*dvp_dzp*dzpdQ
                                              +cn_th*dn_th*(dElKdkth*dkthdzp*dzpdQ*GV_Omega_th + GV_Kcomp_th*dOmegaTHdQ))/PI;

    dVrdr = one_minus_E2*( (GV_R_p_o - GV_R_peri)*(GV_R_p_o-GV_R_4)*(GV_R_apo+GV_R_peri)
                         + (GV_R_apo - GV_R_p_o)*(GV_R_p_o - GV_R_3)*((GV_R_p_o-GV_R_4)+(GV_R_p_o-GV_R_peri)  ) );
    dVthdcth = sin(2.0*theta)*GV_Lz*epsilon_o*(GV_Z_plus+GV_Z_minus-2.0*SQR(cosTh));

/*===========================
    Osculating Conditions
 ============================*/
    
/*--- Derivative of r wrt Wr ----*/
    denom = SQR(sn_r)*(GV_R_apo - GV_R_peri)-(GV_R_apo - GV_R_3);
    
    drdw = 2*sign_r*sn_r*cn_r*dn_r*(GV_R_apo - GV_R_peri)*(GV_R_apo - GV_R_3)*(GV_R_peri - GV_R_3)*GV_Kcomp_r;
    drdw /= SQR(denom)*PI;
    
/*--- Derivative of cos(theta) wrt Wth ----*/
    dcosthdw = 2*sign_th*sqrt(GV_Z_minus)*cn_th*dn_th*GV_Kcomp_th;
    dcosthdw /= PI;
    
   // printf("%d\n",sign_th);
   // printf("%4.6e   %4.6e\n", drdw, dcosthdw);
    
    
    if( (drdw == 0)||(dcosthdw == 0) ){    printf("\n NUL AV DERIVATIVES!\n"); exit(0); }
    
    //printf("%4.6e   %4.6e\n", drdw, dcosthdw);
    
    TOL =1.0E-6;
    
    *dwrdl = drdr3*(dr3dE*GV_E_dot + dr3dQ*GV_Q_dot) + drdr4*(dr4dE*GV_E_dot + dr4dQ*GV_Q_dot);
    *dwrdl /= -drdw;
    
    *dwthdl  = dcosthdzp*(dzpdE*GV_E_dot + dzpdL*GV_Lz_dot + dzpdQ*GV_Q_dot);
    *dwthdl /= -dcosthdw;
    
// (1) First set of oscullating conditions: Singular at ∂r/∂w = 0
// NOTE: On resonace the turnining points for r and theta don't coincide!!
// (2) Second set of oscullating conditions: Singular at ∂V/∂q_µ = 0 NOTE: CHECK THE 2 FACTOR IN FRON OT F

/*  if( (fabs(drdw) > TOL)&&(fabs(dcosthdw) > TOL))
  {
    *dwrdl = drdr3*(dr3dE*GV_E_dot + dr3dQ*GV_Q_dot) + drdr4*(dr4dE*GV_E_dot + dr4dQ*GV_Q_dot);
    *dwrdl /= -drdw;
  
    *dwthdl  = dcosthdzp*(dzpdE*GV_E_dot + dzpdL*GV_Lz_dot + dzpdQ*GV_Q_dot);
    *dwthdl /= -dcosthdw;
  }
  else if(fabs(drdw) > TOL)
  {
    *dwrdl = drdr3*(dr3dE*GV_E_dot + dr3dQ*GV_Q_dot) + drdr4*(dr4dE*GV_E_dot + dr4dQ*GV_Q_dot);
    *dwrdl /= -drdw;
    
    *dwthdl  = GV_Omega_th*2.0*( -GV_pnForce[2] + d2costhdlE*GV_E_dot + d2costhdlL*GV_Lz_dot + d2costhdlQ*GV_Q_dot );
    *dwthdl /= dVthdcth;
  }
  else if(fabs(dcosthdw) > TOL)
  {
    *dwrdl  = GV_Omega_r*2.0*(  -GV_pnForce[1] + d2rdlE*GV_E_dot + d2rdlQ*GV_Q_dot );
    *dwrdl /= dVrdr;

    *dwthdl  = dcosthdzp*(dzpdE*GV_E_dot + dzpdL*GV_Lz_dot + dzpdQ*GV_Q_dot);
    *dwthdl /= -dcosthdw;
   }
   else
   {
    
    *dwrdl  = GV_Omega_r*2.0*(  -GV_pnForce[1] + d2rdlE*GV_E_dot + d2rdlQ*GV_Q_dot );
    *dwrdl /= dVrdr;
    
    *dwthdl  = GV_Omega_th*2.0*( -GV_pnForce[2] + d2costhdlE*GV_E_dot + d2costhdlL*GV_Lz_dot + d2costhdlQ*GV_Q_dot );
    *dwthdl /= dVthdcth;
   }
  */  
//    double  M = 1;//GV_Mass_in_MSun*MASS_SUN*(G_N/GC_c_to_the_3);
    
    *dwrdl /= (2*PI);
    *dwthdl/= (2*PI);
    
/*--- "Geodesic" osculating evolution---*/

    //--Radial geodesic equation
    drdl = drdw*(*dwrdl);
    
    //--Polar geodesic equation
    dcosthdl = dcosthdw*(*dwthdl);
    
      
/*----- This is the end -----*/
    return 0;
}