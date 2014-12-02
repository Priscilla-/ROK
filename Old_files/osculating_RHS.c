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

//int analytic_r_costh(double, double, double *, double *);
void JacobianEllipticFunctions(double, double, double *, double * , double * );
int osculating_turning_points(double *, double *,double *,double *,double *, double *);
int  PN_ELzQdot(double *, double *, double *,double ,double);


int osculating_RHS(double *dWrodt, double *dWTHodt, double *dphidt, double Wr, double Wth)
{
  double drpd[3],drad[3],dr3d[3],dr4d[3],dzmd[3],dzpd[3];// Turning points derivatives
  double drpdE, dradE, drpdLz, dradLz, drpdQ, dradQ;     // "
  double dr3dE, dr4dE, dr3dLz, dr4dLz, dr3dQ, dr4dQ;     // "
  double dzmdE, dzpdE, dzmdLz, dzpdLz, dzmdQ, dzpdQ;     // "
  double drdr3, drdr4,dcosthdzp;            // Derivatives of r wrt the r3 and r4 turning points, and derivative of Costh wrt z+
  double dkrdr3, dkrdr4, dkthdzp;           // Derivatives of the argument of the eliptical function
  double dzpdE,dzpdQ;                       // Derivatives of the polar turning pints
  double one_minus_E2, Sone_minus_E2;       // Auxiliar varialbles
  double dElKdkr,dElKdkth;                  // Derivatives complete elliptical functions
  double r3_root, cos2chi;                  // Aux variables
  double r2pa2,delta,sigma,SQr2pa2;         //   "
  double denom_phi,numer_phi, denom;        //   "
  double epsilon_o,z1;                      // Aux variable
  double q1,q2,q3,q4,q5,q6;
  double varphi_r,varphi_th, dvp_dzp;       // Argument Jacobi's elliptic functions
  double dvp_r3, dvp_r4;                    // Derivatives argument Jacobi elliptical function
  double dVrdr, dVthdcth;                   // Derivatives of the radial and polar potentials
  double sn_r, sn_th;                       // Jacobi's elliptic functions
  double cn_r, cn_th;                       // Jacobi's elliptic functions
  double dn_r, dn_th;                       // Jacobi's elliptic functions
  double drdw, dcosthdw;                    // Radial and plolar derivatives wrt their fundamental frequencies
  double A, B, dAdr3, dAdr4,dBdr3, dBdr4;   // denominator and nominator of radial coordinates and derivatives
  double g1,N,D, dNdr3, dNdr4, dDdE, dDdQ;  // Aux variables
  double dOmegaRdE,dOmegaRdQ;               // Derivatives of the fundamental frequencies wrt constants of the motion
  double dOmegaTHdE, dOmegaTHdL, dOmegaTHdQ;//  "
  double d2rdlE,d2rdlQ;                     // Derivatives of the radial "geodesic eq." wrt constants of the motion
  double d2costhdlE,d2costhdlQ,d2costhdlL;  // Derivatives of the polar "geodesic eq" wrt constants of the motion
  double TOL;                               // Tolerance set to choose the oscullating condition to be applied.
  double cos_psi, cos_chi, rp, costhp;
  int sign_r, sign_th;                      // Aux variables
  
  cos_psi = cos(PSI_NEW);
  cos_chi = COSTH_p_o;
  
  
/*----- Computing dphi_dt -----*/  
  q1 = 1.0 + ECCENTRICITY*cos_psi;
  q2 = P_p/q1;                     // r = p M / ( 1 + e cos(psi) )
  
  q3 = P_p - P3_p - ECCENTRICITY*(P_p + P3_p*cos_psi);  // p - p3 - e ( p + p3 cos(psi) )
  q4 = P_p - P4_p + ECCENTRICITY*(P_p - P4_p*cos_psi);  // p - p4 + e ( p - p4 cos(psi) )
  
  q5 = SQR(q2) + SQR(SPIN);       // r^2 + a^2
  delta = q5 - 2.0*q2;            // Delta = r^2 - 2 M r + a^2
  
  q6 = ENERGY*(SQR(q5)/delta-SQR(SPIN)) - 2.0*SPIN*q2*LZ/delta;  // gamma = E [ ( r^2 + a^2 )^2 / Delta - a^2 ] - 2 M a r Lz / Delta
  
  z1 = Z_minus*SQR(cos_chi);            // z = z_minus cos(chi)^2 = cos(theta)^2
  
  numer_phi = ( (1.0/(1.0-z1) - SQR(SPIN)/delta)*LZ + 2.0*SPIN*ENERGY*q2/delta );
  denom_phi = q6 + SQR(SPIN)*ENERGY*z1;
  
//---Azimutal geodesic equation
  *dphidt = numer_phi/denom_phi;
  
/*----- Computing variations of the turning points -----*/
  osculating_turning_points(&drpd,&drad,&dr3d,&dr4d,&dzmd,&dzpd);
  
  drpdE = drpd[0]; dradE = drad[0]; drpdLz = drpd[1]; dradLz = drad[1]; drpdQ = drpd[2]; dradQ = drad[2];
  dr3dE = dr3d[0]; dr4dE = dr4d[0]; dr3dLz = dr3d[1]; dr4dLz = dr4d[1]; dr3dQ = dr3d[2]; dr4dQ = dr4d[2];
  dzmdE = dzmd[0]; dzpdE = dzpd[0]; dzmdLz = dzmd[1]; dzpdLz = dzpd[1]; dzmdQ = dzmd[2]; dzpdQ = dzpd[2];
  
//---Computing derivatives of the arguments of the Elliptic functions 
  dkrdr3 =  (R_apo - R_peri)*(R_apo-R_4);
  dkrdr3 /= 2.0*KR*(R_peri-R_4)*SQR(R_apo - R_3);
  
  dkrdr4 = (R_apo - R_peri)*(R_3-R_peri);
  dkrdr4 /= 2.0*KR*(R_apo - R_3)*SQR(R_peri-R_4);
  
  dkthdzp =-Z_minus/(2.0*KTH*SQR(Z_plus));
  
  
//---Computing the partial derivatives of the radial and polar turning points wrt constants of motion
  one_minus_E2 = (1.0-SQR(ENERGY));
  Sone_minus_E2 = SQR(one_minus_E2);
  epsilon_o = SQR(SPIN)*(1.0-SQR(ENERGY))/SQR(LZ);
  
  r3_root = sqrt(SQR(2.0/one_minus_E2 - (R_peri+R_apo)) -4.0*SQR(SPIN)*Q_CONSTANT/(one_minus_E2*R_apo*R_peri));
  
  dr3dE = 2.0*ENERGY*(1.0 + (1.0/one_minus_E2 - SQR(SPIN)*Q_CONSTANT/(R_peri*R_apo))/r3_root  );
  dr3dE /= Sone_minus_E2;
  
  dr4dE = SQR(SPIN)*Q_CONSTANT*(2.0*ENERGY*R_3/Sone_minus_E2-dr3dE/one_minus_E2);
  dr4dE /= R_peri*R_apo*SQR(R_3);
  
  dr3dQ = -SQR(SPIN);
  dr3dQ /= one_minus_E2*R_peri*R_apo*r3_root;
  
  dr4dQ = SQR(SPIN)*(R_3 -Q_CONSTANT*dr3dQ);
  dr4dQ /= one_minus_E2*R_peri*R_apo*SQR(R_3);
  
  dzpdE = 2.0*Q_CONSTANT*ENERGY;
  dzpdE /=Z_minus*SQR(SPIN*(1.-SQR(ENERGY)));
  
  dzpdQ = 1.0/(SQR(SPIN)*Z_minus*(1.-SQR(ENERGY)));
  
  // Sanity checks
  /*double r3= 0.5*(2/one_minus_E2-(R_peri+R_apo) +r3_root);
   double r4  = SQR(SPIN)*Q_CONSTANT/(one_minus_E2*R_apo*R_peri*r3);
   double r3_p = P3_p/(1-ECCENTRICITY);
   double r4_p = P4_p/(1+ECCENTRICITY);
   printf("r3=%4.6e, r3_p=%4.6e, R3=%4.6e, r4= %4.6e, r4_p= %4.6e, R4= %4.6e\n", r3, r3_p, R_3, r4, r4_p, R_4);
   */
  
//---- Derivatives of the Complete elliptical integrals of first kind K(k_r), K(k_th) 
  dElKdkr  = Ecom_r/(KR*(  1.0-SQR(KR))) - Kcom_r/KR;
  dElKdkth = Ecom_th/(KTH*(1.0-SQR(KTH)))- Kcom_th/KTH;
  
//---- Derivatives of the fundamental frequencies  
  dOmegaRdE = -SQR(PI)*(2.0*ENERGY*(R_apo-R_3)*(R_peri-R_4) +(1.0-SQR(ENERGY))*( (R_peri-R_4)*dr3dE + (R_apo-R_3)*dr4dE));
  dOmegaRdE /= OMEGA_r*8.0*SQR(Kcom_r);
  dOmegaRdE += -OMEGA_r*dElKdkr*(dkrdr3*dr3dE + dkrdr4*dr4dE)/Kcom_r;
  
  dOmegaRdQ = -SQR(PI)*( (1.-SQR(ENERGY))*( (R_peri-R_4)*dr3dQ + (R_apo-R_3)*dr4dQ));
  dOmegaRdQ /= OMEGA_r*8.0*SQR(Kcom_r);
  dOmegaRdQ += -OMEGA_r*dElKdkr*(dkrdr3*dr3dQ + dkrdr4*dr4dQ)/Kcom_r;
  
  dOmegaTHdE =  PI*LZ*(-Z_plus*(2.0*SPIN2*ENERGY)/SQR(LZ) + epsilon_o*dzpdE );
  dOmegaTHdE /= 4.0*Kcom_th*sqrt(epsilon_o*Z_plus);
  dOmegaTHdE += -OMEGA_th*dElKdkth*dkthdzp*dzpdE/Kcom_th;
  
  dOmegaTHdQ =  PI*LZ*epsilon_o*dzpdQ;
  dOmegaTHdQ /= 4.0*Kcom_th*sqrt(epsilon_o*Z_plus);
  dOmegaTHdQ += -OMEGA_th*dElKdkth*dkthdzp*dzpdQ/Kcom_th;
  
  dOmegaTHdL =  PI*LZ*(-Z_plus*(2.0*SPIN2*ENERGY)/SQR(LZ));
  dOmegaTHdL /= 4.0*Kcom_th*sqrt(epsilon_o*Z_plus);
  dOmegaTHdL += (0.5*PI*sqrt(epsilon_o*Z_plus))/Kcom_th;
  
/*--- Arguments of the Jacobi elliptic functions sn(varphi,k), cn(varphi,k), dn(varphi,k) ----*/
  sign_r =0;
  sign_th=0;
    
    // printf(" Wr=%4.6e  Wth =%4.6e\n",Wr ,Wth);
    
// Radial argument
  if( (Wr>= 0)&&(Wr<=PI) )
  {
      varphi_r = Wr*Kcom_r/PI;
      sign_r = 1;
  }
  else if( (Wr>= PI)&&(Wr<=2.*PI) )
  {
      varphi_r = (2.*PI-Wr)*Kcom_r/PI;
      sign_r = -1;
  }
  else
      printf("Wr out of bounds!!\n");
  
  // Polar argument
  if( (Wth>= 0)&&(Wth<=0.5*PI) )
  {
      varphi_th = Wth*0.5*Kcom_th;
      sign_th = 1;
  }
  else if( (Wth>= 0.5*PI)&&(Wth<=1.5*PI) )
  {
      varphi_th = (PI-Wth)*2.*Kcom_th/PI;
      sign_th =-1;
  }
  else if( (Wth>= 0.5*PI)&&(Wth<=PI) )
  {
      varphi_th = (Wth-2.*PI)*2.*Kcom_th/PI;
      sign_th = 1;
  }
  else
      printf("Wth out of bounds!!\n");
  
  
/*---- Derivatives of the arguments of the Jacobi elliptical function sn(varphi,k), cn(varphi,k), dn(varphi,k) ---*/
  // Radial argument
  if( (Wr>= 0.)&&(Wr<=PI) )
  {
      dvp_r3 = Wr*dElKdkr*dkrdr3;
      dvp_r3 /= PI;
      
      dvp_r4 = Wr*dElKdkr*dkrdr4;
      dvp_r4 /= PI;
  }
  else if( (Wr>= PI)&&(Wr<=2.*PI) )
  {
      dvp_r3 = (2.0*PI-Wr)*dElKdkr*dkrdr3;
      dvp_r3 /= PI;
      
      dvp_r4 = (2.0*PI-Wr)*dElKdkr*dkrdr4;
      dvp_r4 /= PI;
  }
  else
    printf("dvp_r out of bounds!!\n");
  
  // Polar argument
  if( (Wth>= 0.)&&(Wth<=0.5*PI) )
  {
    dvp_dzp = Wth*2.0*dElKdkth*dkthdzp/PI;
    // printf("(1)\n");
  
  }
  else if( (Wth>= 0.5*PI)&&(Wth<=1.5*PI) ){
      dvp_dzp = (PI-Wth)*2.0*dElKdkth*dkthdzp/PI;
    // printf("(2)\n");
  }
  else if( (Wth>= 1.5*PI)&&(Wth<=PI) ){
      dvp_dzp = (Wth-2.*PI)*2.*dElKdkth*dkthdzp/PI;
    //  printf("(3)\n");
  }
  else
    printf("dvp_th out of bounds!!\n");

/*---Jacobi elliptic functions----*/
 JacobianEllipticFunctions(varphi_r, KR,  &sn_r,  &cn_r ,  &dn_r);
 JacobianEllipticFunctions(varphi_th,KTH, &sn_th, &cn_th , &dn_th);

/*---- Partial derivatives of r(w_b)  [r = A/B ] wrt the turning points r3, r4----*/
  A = R_3*(R_apo - R_peri)*SQR(sn_r) - R_peri*(R_apo - R_3);
  B = (R_apo - R_peri)*SQR(sn_r) -(R_apo - R_3);
  
  dAdr3 = (R_apo - R_peri)*(SQR(sn_r) + R_3*2.0*sn_r*cn_r*dn_r*dvp_r3) + R_peri;
  dAdr4 = R_3*(R_apo - R_peri)*2.0*sn_r*cn_r*dn_r*dvp_r4;
  
  dBdr3 = (R_apo - R_peri)*2.0*sn_r*cn_r*dn_r*dvp_r3 + 1.0;
  dBdr4 = (R_apo - R_peri)*2.0*sn_r*cn_r*dn_r*dvp_r4;
  
  drdr3 = B*dAdr3-A*dBdr3;
  drdr3 /=SQR(B);
  
  drdr4 = B*dAdr4-A*dBdr4;
  drdr4 /=SQR(B);
  
/*---- Partial derivatives of cos(theta)(Wth) wrt turning points----*/
  dcosthdzp = sqrt(Z_minus)*sqrt(1.-SQR(sn_th))*dn_th*dvp_dzp*(dzpdE + dzpdQ);
  
   // printf("dvp_dzp=%4.6e\n",dElKdkth*dkthdzp/PI);
    
/*---- Computing the second derivatives of the geodesic equations [dr/dλ = N*D/B^2, dcos(theta)[Wth]/dλ] wrt the constant of motion---*/
  
  g1 = (R_peri-R_3)*(R_apo-R_3)*(R_apo-R_peri);
  N  = sn_r*cn_r*dn_r*g1;
  D  = Kcom_r*OMEGA_r/PI;
  
  dNdr3 = (SQR(cn_r*dn_r) - SQR(sn_r)*( SQR(dn_r) + SQR(KR*cn_r) ) )*g1*dvp_r3 - sn_r*cn_r*dn_r*(R_apo+R_peri)*(R_apo-R_peri);
  dNdr4 = (SQR(cn_r*dn_r) - SQR(sn_r)*( SQR(dn_r) + SQR(KR*cn_r) ) )*g1*dvp_r4;
  
  dOmegaRdE = -SQR(PI)*(2.0*ENERGY*(R_apo-R_3)*(R_peri-R_4) +(1.0-SQR(ENERGY))*( (R_peri-R_4)*dr3dE + (R_apo-R_3)*dr4dE));
  dOmegaRdE /= OMEGA_r*8.0*SQR(Kcom_r);
  dOmegaRdE += -OMEGA_r*dElKdkr*(dkrdr3*dr3dE + dkrdr4*dr4dE)/Kcom_r;
  
  dOmegaRdQ = -SQR(PI)*( (1.-SQR(ENERGY))*( (R_peri-R_4)*dr3dQ + (R_apo-R_3)*dr4dQ));
  dOmegaRdQ /= OMEGA_r*8.0*SQR(Kcom_r);
  dOmegaRdQ += -OMEGA_r*dElKdkr*(dkrdr3*dr3dQ + dkrdr4*dr4dQ)/Kcom_r;
  
  dOmegaTHdE =  PI*LZ*(-Z_plus*(2.0*SPIN2*ENERGY)/SQR(LZ) + epsilon_o*dzpdE );
  dOmegaTHdE /= 4.0*Kcom_th*sqrt(epsilon_o*Z_plus);
  dOmegaTHdE += -OMEGA_th*dElKdkth*dkthdzp*dzpdE/Kcom_th;
  
  dOmegaTHdQ =  PI*LZ*epsilon_o*dzpdQ;
  dOmegaTHdQ /= 4.0*Kcom_th*sqrt(epsilon_o*Z_plus);
  dOmegaTHdQ += -OMEGA_th*dElKdkth*dkthdzp*dzpdQ/Kcom_th;
  
  dOmegaTHdL =  PI*LZ*(-Z_plus*(2.0*SPIN2*ENERGY)/SQR(LZ));
  dOmegaTHdL /= 4.0*Kcom_th*sqrt(epsilon_o*Z_plus);
  dOmegaTHdL += (0.5*PI*sqrt(epsilon_o*Z_plus))/Kcom_th;
  
  dDdE = dElKdkr*(dkrdr3*dr3dE + dkrdr4*dr4dE)*OMEGA_r/PI + Ecom_r*dOmegaRdE;
  dDdQ = dElKdkr*(dkrdr3*dr3dQ + dkrdr4*dr4dQ)*OMEGA_r/PI + Ecom_r*dOmegaRdQ;
  
  
  d2rdlE  = (2.0*sign_r*( (dNdr3*dr3dE + dNdr4*dr4dE)*SQR(B) - 2.0*N*B*(dBdr3*dr3dE + dBdr4*dr4dE))*D)/SQR(B*B);
  d2rdlE += N*(dElKdkr*(dkrdr3*dr3dE + dkrdr4*dr4dE )*OMEGA_r + Kcom_r*dOmegaRdE)/(PI*SQR(B));
  
  d2rdlQ  = (2.0*sign_r*( (dNdr3*dr3dQ + dNdr4*dr4dQ)*SQR(B) - 2.0*N*B*(dBdr3*dr3dQ + dBdr4*dr4dQ))*D)/SQR(B*B);
  d2rdlQ += N*(dElKdkr*(dkrdr3*dr3dQ + dkrdr4*dr4dQ )*OMEGA_r + Kcom_r*dOmegaRdQ)/(PI*SQR(B));
  
  
  d2costhdlE = 2.0*sign_th*sqrt(Z_minus)*( (-(sn_th*SQR(dn_r) + SQR(KTH*cn_th)*sn_th)*Kcom_th*OMEGA_th*dvp_dzp*dzpdE
                                 +cn_th*dn_th*(dElKdkth*dkthdzp*dzpdE*OMEGA_th + Kcom_th*dOmegaTHdE)))/PI;
  d2costhdlL = 2.0*sign_th*sqrt(Z_minus)*(-(sn_th*SQR(dn_r)+ Kcom_th*dOmegaTHdL))/PI;
  d2costhdlQ = 2.0*sign_th*sqrt(Z_minus)*(-(sn_th*SQR(dn_r) + SQR(KTH*cn_th)*sn_th)*Kcom_th*OMEGA_th*dvp_dzp*dzpdQ
                               +cn_th*dn_th*(dElKdkth*dkthdzp*dzpdQ*OMEGA_th + Kcom_th*dOmegaTHdQ))/PI;
  
  dVrdr = one_minus_E2*( (R_p_o - R_peri)*(R_p_o-R_4)*(R_apo+R_peri)
                  + (R_apo - R_p_o)*(R_p_o - R_3)*((R_p_o-R_4)+(R_p_o-R_peri)  ) );
  dVthdcth = sin(2.0*acos(COSTH_p_o))*LZ*epsilon_o*(Z_plus+Z_minus-2.0*SQR(COSTH_p_o));
  
/*===========================
   Osculating Conditions
============================*/
  // (1) First set of oscullating conditions: Singular at ∂r/∂w = 0
  // NOTE: On resonace the turnining points for r and theta don't coincide!!
  // (2) Second set of oscullating conditions: Singular at ∂V/∂q_µ = 0 NOTE: CHECK THE 2 FACTOR
  
//--- Derivative of r wrt Wr ----
   denom = SQR(sn_r)*(R_apo - R_peri)-(R_apo - R_3);
   
   drdw = 2.*sign_r*sn_r*cn_r*dn_r*(R_apo - R_peri)*(R_apo - R_3)*(R_peri - R_3)*Kcom_r;
   drdw /= SQR(denom)*PI;
   
   //--- Derivative of cos(theta) wrt Wth ---
   dcosthdw = 2.*sign_th*sqrt(Z_minus)*cn_th*dn_th*Kcom_th;
   dcosthdw /= PI;
   
   TOL =1.0E-6;
   
   sigma = SQR(R_p_o) + SPIN2*SQR(COSTH_p_o);
   
   if (fabs(drdw) > TOL)
   {
   *dWrodt = drdr3*(dr3dE*ENERGY_dot + dr3dQ*Q_CONSTANT_dot) + drdr4*(dr4dE*ENERGY_dot + dr4dQ*Q_CONSTANT_dot);
   *dWrodt /= -drdw;
   }
   else
   {
   *dWrodt  = OMEGA_r*( -PNFcomp[1] + sigma*2.0*(d2rdlE*ENERGY_dot + d2rdlQ*Q_CONSTANT_dot) );
   *dWrodt /= dVrdr;
   }
   
   if (fabs(dcosthdw) > TOL)
   {
   *dWTHodt  = dcosthdzp*(dzpdE*ENERGY_dot + dzpdQ*Q_CONSTANT_dot);
   *dWTHodt /= -dcosthdw;
   }
   else
   {
   *dWTHodt  = OMEGA_th*( -PNFcomp[2] + sigma*2.0*(d2costhdlE*ENERGY_dot + d2costhdlL*LZ_dot + d2costhdlQ*Q_CONSTANT_dot) );
   *dWTHodt /= dVthdcth;
   }
    
/*----- This is the end -----*/
  return 0;
}


// Alternative Forumulations

/*
 double dpsi_o_dt, dchi_o_dt;
 double dwrdpsi, dwthdchi;
 double dwrdr, dwthdcosth;
 double J_psi,denom;
 //double dpsidr, dchidr;
 
 J_psi = (1.-SQR(ENERGY))*(1.-SQR(ECCENTRICITY)) + 2.*(1. - SQR(ENERGY) - (1. -SQR(ECCENTRICITY)/P_p) )*(1.+ECCENTRICITY*cos(psi))
 +( (1.-SQR(ENERGY))*(3.+SQR(ECCENTRICITY))/(1.-SQR(ECCENTRICITY)) -4./P_p
 +(1.-SQR(ECCENTRICITY))/SQR(P_p)*( SPIN2*(1.-SQR(ENERGY)) + SQR(LZ) + Q_CONSTANT))*SQR(1. +ECCENTRICITY*cos(psi));
 
 dwrdpsi = (1.-SQR(ECCENTRICITY))*OMEGA_r;
 dwrdpsi /= P_p*sqrt(J_psi);
 
 denom = SQR(sn_r)*(R_apo - R_peri)-(R_apo - R_3);
 
 dwrdr  = SQR(denom)*PI;
 dwrdr /= 2.*sign_r*sn_r*cn_r*dn_r*(R_apo - R_peri)*(R_apo - R_3)*(R_peri - R_3)*Kcom_r;
 
 *dpsi_o_dt  = dwrdr*( drdr3*(dr3dE*ENERGY_dot + dr3dQ*Q_CONSTANT_dot) + drdr4*(dr4dE*ENERGY_dot + dr4dQ*Q_CONSTANT_dot));
 *dpsi_o_dt /= dwrdpsi;
 
 
 dwthdchi = PI;
 dwthdchi /= 2.Kcom_th*(1.-SQR(KTH*cos(chi)));
 
 dwthdcosth = PI;
 dwthdcosth /=2.*sign_th*sqrt(Z_minus)*cn_th*dn_th*Kcom_th;
 
 dcosthdw = 2.*sign_th*sqrt(Z_minus)*cn_th*dn_th*Kcom_th;
 dcosthdw /= PI;
 
 *dchi_o_dt = dwthdcosth*dcosthdzp*(dzpdE*ENERGY_dot + dzpdQ*Q_CONSTANT_dot);
 *dchi_o_dt/= dwthdchi;
 
 //  printf("dcosthdzp=%4.6e  OMEGA_th =%4.6e  PNFcomp[2]=%4.6e  Wth =%4.6e\n",dcosthdzp,OMEGA_th, PNFcomp[2] ,Wth);
 
 */


// (1) First set of oscullating conditions: Singular at ∂r/∂w = 0
// NOTE: On resonace the turnining points for r and theta don't coincide!!
// (2) Second set of oscullating conditions: Singular at ∂V/∂q_µ = 0 NOTE: CHECK THE 2 FACTOR

/*  //--- Derivative of r wrt Wr ----
 denom = SQR(sn_r)*(R_apo - R_peri)-(R_apo - R_3);
 
 drdw = 2.*sign_r*sn_r*cn_r*dn_r*(R_apo - R_peri)*(R_apo - R_3)*(R_peri - R_3)*Kcom_r;
 drdw /= SQR(denom)*PI;
 
 //--- Derivative of cos(theta) wrt Wth ---
 dcosthdw = 2.*sign_th*sqrt(Z_minus)*cn_th*dn_th*Kcom_th;
 dcosthdw /= PI;
 
 TOL =1.0E-6;
 
 sigma = SQR(R_p_o) + SPIN2*SQR(COSTH_p_o);
 
 if (fabs(drdw) > TOL)
 {
 *dWrodt = drdr3*(dr3dE*ENERGY_dot + dr3dQ*Q_CONSTANT_dot) + drdr4*(dr4dE*ENERGY_dot + dr4dQ*Q_CONSTANT_dot);
 *dWrodt /= -drdw;
 }
 else
 {
 *dWrodt  = OMEGA_r*( -PNFcomp[1] + sigma*2.0*(d2rdlE*ENERGY_dot + d2rdlQ*Q_CONSTANT_dot) );
 *dWrodt /= dVrdr;
 }
 
 if (fabs(dcosthdw) > TOL)
 {
 *dWTHodt  = dcosthdzp*(dzpdE*ENERGY_dot + dzpdL*LZ_dot + dzpdQ*Q_CONSTANT_dot);
 *dWTHodt /= -dcosthdw;
 }
 else
 {
 *dWTHodt  = OMEGA_th*( -PNFcomp[2] + sigma*2.0*(d2costhdlE*ENERGY_dot + d2costhdlL*LZ_dot + d2costhdlQ*Q_CONSTANT_dot) );
 *dWTHodt /= dVthdcth;
 }
 */

// ALternative formulation check

/*
 double dpsidwr, dchidwth;
 double dpsidr, dchidcosth;
 double J;  //Schmidt’s J function
 double arg, costh;
 
 J = (1.-SQR(ENERGY))*(1.-SQR(ECCENTRICITY)) + 2.*(1. - SQR(ENERGY) - (1. -SQR(ECCENTRICITY)/P_p) )*(1.+ECCENTRICITY*cos(PSI_NEW))
 +( (1.-SQR(ENERGY))*(3.+SQR(ECCENTRICITY))/(1.-SQR(ECCENTRICITY)) -4./P_p
 +(1.-SQR(ECCENTRICITY))/SQR(P_p)*( SPIN2*(1.-SQR(ENERGY)) + SQR(LZ) + Q_CONSTANT))*SQR(1. +ECCENTRICITY*cos(PSI_NEW));
 
 dpsidwr = P_p*sqrt(J);
 dpsidwr /=(1.-SQR(ECCENTRICITY))*OMEGA_r;
 
 arg = (P_P/R_p_o - 1.)/ECCENTRICITY;
 dpsidr = P_P/((1.-SQR(arg))*(ECCENTRICITY*SQR(R_p_o)));
 
 *dWrodt = dpsidr*( drdr3*(dr3dE*ENERGY_dot + dr3dQ*Q_CONSTANT_dot) + drdr4*(dr4dE*ENERGY_dot + dr4dQ*Q_CONSTANT_dot));
 *dWrodt/= dpsidwr;
 
 dchidwth = 2.Kcom_th*(1.-SQR(KTH*cos(CHI_NEW)));
 dchidwth /= PI;
 
 costh = qrt(Z_minus)*sn_th;
 
 dchidcosth = -1./(1.-SQR(costh/Z_minus))*Z_minus;
 
 *dWTHodt = dchidcosth*dcosthdzp*(dzpdE*ENERGY_dot + dzpdL*LZ_dot + dzpdQ*Q_CONSTANT_dot);
 *dWTHodt /=dchidwth;
 
 */










