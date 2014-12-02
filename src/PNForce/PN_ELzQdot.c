/*==============================================================================
 Determines contravariant, f[], and covariant, fcov[], self-force at a point
 where the velocity is u[] and the covariant velocity ucov[].
 The vecor params contains various necessary variables that are computed within
 the Kerrderivs code. These are, in order
 r, cos(theta), sin(theta), E, Lz, Q, cosiota, z, delta, sigma, chi
 -------------------------------------------------------------------------------
 PN Force [only depends on the two orbital phases psi and chi] facilitated by 
 Eanna Flanagan.
 
 Here we follow the notation and formulation of Gair et al: PRD 83, 044037 [1]
 
 NOTE: Most of the comments here belong to the Eanna & Tanja's routine
 -------------------------------------------------------------------------------
 Last Update : 23.05.13
 ============================================================================*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
 
#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"

#include "physical_quantities.h"
#include "macros_2.h"

//int dchi_dpsi_dphi_dt(double *, double *, double, double);

int PN_ELzQdot(double *dEndt, double *dLzdt, double *dQdt, double *dxdt, double *x, double r_BL, double cosTh_BL,double dpsidt, double dchidt,double dphidt, double dtaudt)
{
  //double dChidt, dPhidt, dPsidt,dRdt, dTaudt, dThetadt;   // derivatives of the BL coordinates in terms of the BL time t
  double psi, chi,tau, p, ecc,iota,M;
  double ac1,U0,U1,U2,Ut,bara2p,ac2,Vt,an,Ra,Ia;
  double E,Lz2,K,cosiota;           // Energy, azimuthal angular momentum,  Carter constant, cos inclination angle [ in Jon's code E=params[3],
  double zeta;                       // z = z_minus cos(chi)^2 = cos(theta)^2 [in Jon's code z=params[7] ]
  double delta,sigma;                // Auxiliar variables [ in Jon's code delta=params[8], sigma=params[9] ]
  double Ffnc, calH;                 // Ffnc = eq. 3.32 [1], calH = eq. 3.31 [1]
  double unEF,RuEF,IuEF,ulEF;        // Four velocity components in terms of the Kinnersley tetrad
  double u[4], ucov[4];              // Contravariant u^{\mu} and Covariant u_{\mu} velocity components
  double f[4], fcov[4];              // Contravariant u^{\mu} and Covariant f_{\mu} force components
  double gtt, grr, gthth, gphiphi, gtphi;// Kerr metric components g_{\mu}{\nu}
  double detg;                       // Kerr metric determinant
  // PN SF aux variables:
  double rEF, rsq, rth;              // rEF= r + a^2/(4*r), rEF^2, rEF^3
  double r4, r5, r6, r7, r8, r9;     // rEF^4,...,rEF^9
  double cottheta, cosectheta, cos2iota;
  double cos2theta, j1, j2;
  double drdt, dthetadt;
  double dtilderdtildet;
  double dthetadtildet;
  double EN, EN2, Kqp2;
  double al;
  double ac10, ac11, ac12;
  double a2p0, a2p1, a2p2;
  double rootdelta, rootsigma, baru1;
  double baru2,chi2,chiEF;
  double bara1p, bara1, bara2, alll;
  double r2, r2pa2,lhs1,lhs2;
  double coshpsi, sinhpsi, cos_theta, sin_theta,sin_psi,sin_chi, cos_psi,cos_chi;
  double gamma, Kqp,opecp;
  double dtau_dpsi, dtau_dchi, dr_dpsi, dth_dchi;
  double drdtau, dcosthdtau;
    int i;

  EPSILON = MAG_SF;
  
  p   = x[3];
  ecc = x[4];
  iota= x[5];
    
/*----- Auxiliary quantities -----*/

  sin_theta = sqrt(1. -SQR(cosTh_BL));
    
  cos_psi = (x[3]/r_BL - 1.)/x[4];
  sin_psi = sqrt(1. - cos_psi*cos_psi);
  
  cos_chi = cosTh_BL/sqrt(Z_minus);
  sin_chi = sqrt(1. -cos_chi*cos_chi);
  opecp = 1. + ecc*cos_psi;
  
  dr_dpsi  = p*ecc*sin_psi/(opecp*opecp);
  dth_dchi = sqrt(fabs(Z_minus))*sin_chi/sin_theta;
   
/*=== Identifying variables with Eanna and Tanja's ones  ===*/ // WARNING: in Fujita's notation Q --> C
/*---- Contravariant Four velocity components (dX/dTau)----*/
  M = MASS_in_MSUN*MASS_SUN*(GN/CSPEED3);


  u[0] = 1./dtaudt;           //dt/dtau
  u[1] = dr_dpsi*dpsidt*u[0]; //dr/dtau
  u[2] = dth_dchi*dchidt*u[0];//dtheta/dtau
  u[3] = dphidt*u[0];         //dphi/dtau
   
// Need this for the EF's PN-SF
  drdt =  u[1]/u[0];
  dthetadt = u[2]/u[0];

 // printf(" drdt =%4.6e, dthetadt =%4.6e \n\n\n",drdt, dthetadt );
      
/*---- Identifying Global variables with Eanna and Tanja's ones  ----*/ // WARNING: in Fujita's notation Q --> C
  cos_theta = cosTh_BL;
    
  zeta    = SQR(cos_theta);
  r2      = r_BL*r_BL;
  r2pa2   = r2 + SPIN*SPIN;
  
  sigma = r2  + SQR(SPIN*cos_theta);
  delta = r2  + SPIN2 -2.0*r_BL;
  gamma = ENERGY*(SQR(r2pa2)/delta-SPIN2) - 2.0*SPIN*r_BL*LZ/delta;    // gamma = E [ ( r^2 + a^2 )^2 / Delta - a^2 ] - 2 M a r Lz / Delta

  E   = ENERGY;
  EN  = E - 1.0;
  EN2 = EN*EN;

  Lz2 = LZ*LZ;
  
  K   = C_CONSTANT;
  Kqp = K + 2.0*SPIN*LZ*E + SPIN2*(-E*E + (E-1.0));//k +2.0*a*LZ*E - a*a*E*E + a*a*(E-1.0);
  Kqp2= SQR(Kqp);
    
  cosiota = cos(iota);
    
  Ffnc = (SQR(r_BL) + SPIN2 )*E - SPIN*LZ;
  calH =  LZ - SPIN*(1.0 - zeta )*E;
    
/*---- EF's PN-SF auxiliar quantities ----*/
  rEF = r_BL + SPIN2/(4.*r_BL);
  rsq = SQR(rEF);
  rth = rsq*rEF;
  r6  = rth*rth;
  r4  = rsq*rsq;
  r5  = rEF*r4;
  r7  = rth*r4;
  r8  = r4*r4;
  r9  = r4*r5;
    
  cottheta   = cos_theta/sin_theta; // cth/sth;
  cosectheta = 1.0/sin_theta;       // 1.0/sth;
  cos2iota   = 2.0*cosiota*cosiota-1.0;
  cos2theta  = 2.0*SQR(cos_theta)-1.0;  //  2.0*cth*cth-1.0;
    
  j1 = SPIN2/rsq;
  j2 = j1*cos2theta;
    
/*----- Computing the  Boyer-Lindquist radial and polar Velocities wrt the Boyer-Lindquist time T_p -----*/
 /* if ( (chi_2pi >= 0.0) && (chi_2pi < PI) )
    dThetadt = sqrt( (Z_minus-zeta)/(1.0-zeta) )*dChidt;// /M;
  else if ( (chi_2pi >= PI) && (chi_2pi < 2.0*PI) )
    dThetadt = - sqrt( (Z_minus-zeta)/(1.0-zeta) )*dChidt;// /M;

  dRdt= (r2/pp)*ecc*sin(psi)*dPsidt;
   */ 

  dtilderdtildet = (1.0 - j1/4.0) * (1.0+j2/2.0) * drdt;
  dthetadtildet = (1.0+j2/2.0) * dthetadt;
  
/*---- Metric components ----*/
  gtt = -1.0 + 2.0*r_BL/sigma;
  gtphi = -2.0*SPIN*SQR(sin_theta)*r_BL/sigma;
  grr   = sigma/delta;
  gphiphi= ( SQR(r2pa2)-delta*SPIN2*SQR(sin_theta) )*SQR(sin_theta)/sigma;
  gthth  = sigma;
  detg = SQR(gtphi)- gtt*gphiphi;
  
 /*---- Covariant Four velocity components----*/
  ucov[0] = gtt*u[0] + gtphi*u[3];
  ucov[1] = grr*u[1];
  ucov[2] = gthth*u[2];
  ucov[3]=  gphiphi*u[3] + gtphi*u[0];

/**---- Four velocity components in therm of the Kinnersley tetrad ----**/
  ulEF = ucov[1]-Ffnc/delta;
  unEF = -(Ffnc + ucov[1]*delta)/(2.*sigma);
  RuEF = r_BL*ucov[2]/sigma + SPIN*cos_theta*calH/(sin_theta*sigma);
  IuEF = SPIN*cos_theta*ucov[2]/sigma - r_BL*calH/(sigma*sin_theta);
    
/*---- Computation of the PN self-force ----*/
  ac10 = (16.0*(15.0*Kqp + 2.0*rEF*(5.0 + 3.0*rEF*EN)))/(15.0*r5);
  ac10 *= dtilderdtildet;
  ac11 = (-8.0*LZ*(105.0*Kqp + 8.0*rEF*(-10.0 + 3.0*rEF*EN)))/(15.0*r7);
  ac11 *= dtilderdtildet;
  ac12 = - (174.0*Kqp*cos_theta*sin_theta*dthetadtildet/r6) + (28.0*cos_theta*sin_theta*dthetadtildet)/(r5)
            + dtilderdtildet*( (128.0*Lz2)/(5.0*r8)
                                - (20.0*Kqp)/r7
                                - (160.0*Lz2)/r7
                                - 58.0/(3.0*r6)
                                - (8.0*EN)/r5
                                - (268.0*Kqp*cos2theta)/r7
                                + (54.0*cos2theta)/r6
                                + (328.0*EN*cos2theta)/(5.0*r5)
                                + (288.0*dtilderdtildet*cos_theta*sin_theta*dthetadtildet)/r4);
  ac1 = ac10 + SPIN*ac11 + SPIN2*ac12;
    
 /* Ut is a^\hat{\theta} times a factor (\sqrt{Sigma}?). */
  U0 = (8.0*(-15.0*Kqp + 2.0*rEF*(10.0 + 9.0*rEF*EN))*dthetadtildet)/(5.0*rth);
  U1 = (-16.0*LZ*(15.0*dtilderdtildet*(-13.0*Kqp + 8.0*rEF*(1.0 + rEF*EN))*cottheta + rEF*(-225.0*Kqp + rEF*(461.0 + 342.0*rEF*EN))*dthetadtildet))/(15.0*r6);
  U2 = (15.0*dtilderdtildet*(-209.0*Kqp + 348.0*Lz2 + 6.0*rEF + 16.0*rsq*EN + (209.0*Kqp - 2.0*rEF*(3.0 + 8.0*rEF*EN))*cos2theta)*      cottheta +rEF*(405.0*Kqp + 1080.0*Lz2 - 832.0*rEF - 648.0*rsq*EN + (4095.0*Kqp - 2.0*rEF*(2815.0 + 2664.0*rEF*EN))*cos2theta)*dthetadtildet)/(15.0*r6);
    Ut = U0 + SPIN*U1 + SPIN2*U2;
    
    /* bara2p is a_perpendicular */
  a2p0 = (8.0*LZ*(-15.0*Kqp + 2.0*rEF*(10.0 + 9.0*rEF*EN))*cosectheta)/(5.0*r6);
  a2p1 = (2.0*cosectheta*(-375.0*Kqp2 + 1800.0*Kqp*Lz2 + 93.0*Kqp*rEF - 3688.0*Lz2*rEF + 790.0*rsq + 420.0*Kqp*rsq*EN - 2736.0*Lz2*rsq*     EN + 936.0*rth*EN + 192.0*r4*EN2 + (375.0*Kqp2 - 3.0*Kqp*rEF*(31.0 + 140.0*rEF*EN) - 2.0*rsq*(395.0 + 468.0*rEF*EN + 96.0*rsq*EN2))*cos2theta + 60.0*rth*dtilderdtildet*(-13.0*Kqp + 8.0*rEF*(1.0 + rEF*EN))*2.*cos_theta*sin_theta*dthetadtildet))/(15.0*r8);
  a2p2 = (LZ*cosectheta*(345.0*Kqp + 360.0*Lz2 - 593.0*rEF - 516.0*rsq*EN + (1185.0*Kqp - rEF*(1601.0 + 1512.0*rEF*EN))*cos2theta - 870.0*rth*dtilderdtildet*2.*cos_theta*sin_theta*dthetadtildet))/(5.0*r8);
  bara2p = a2p0 + SPIN*a2p1 + SPIN2*a2p2;
    
  if(FLAG_incConsPiece == 'y')
  {
    ac1 += Kqp/(2.0*r4) - 5.0/rth - 7.0*EN/rsq;
     //ac2 += dtilderdtildet*( 3.0*Kqp/(2.0*r4) + 5.0/r3 + 7.0*EN/r2);
    Ut += -2.0*dtilderdtildet*dthetadtildet;
     //Vt += -2.0*Lz*dtilderdtildet*cosectheta/r2;
    bara2p += -2.0*LZ*dtilderdtildet*cosectheta/rth;
            
      // now add leading spin terms
    ac1   += 3.0*SPIN*LZ/r4;
    Ut    += -3.0*SPIN*LZ*cottheta/rth;
    bara2p+= 3.0*SPIN*cos_theta*dthetadtildet/rsq;
  }
        
  rootdelta = sqrt(delta);
  rootsigma = sqrt(sigma);
  baru1 = Ffnc/(rootdelta*rootsigma);
    
   // printf("%4.6e %4.6e %4.6e\n",rootdelta,rootsigma,Ffnc);
    //CHECK(baru1<0,"Past directed velocity in computeAcc");
  if(baru1<0) printf("Past directed velocity in computeAcc\n");
    
  baru2 = calH/(rootsigma*sin_theta);   // calH / (rootsigma* sth);
  chi2 = baru1*baru1 - baru2*baru2;
    
    
   // CHECK(chi2 < 0,"Imaginary boost in computeAcc");
  if(chi2 < 0) printf("Imaginary boost in computeAcc\n");
    
  chiEF = sqrt(chi2);

  coshpsi = baru1 / chiEF;
  sinhpsi = baru2 / chiEF;
        
  bara1p = (delta * ac1 * ucov[1] + Ut*ucov[2])/(sigma*chiEF);
  bara1 = coshpsi *bara1p + sinhpsi *bara2p;
  bara2 = sinhpsi *bara1p + coshpsi *bara2p;
  ac2   = - rootsigma*bara1/rootdelta;
  Vt    = rootsigma*bara2;
    
  alll = sigma/delta;
  an = (ac2 - ac1)/(2.0*alll);
  Ra = (r_BL*Ut + SPIN*cos_theta*Vt)/sigma;
  Ia = (SPIN*cos_theta*Ut - r_BL*Vt)/sigma;
        
  an *= EPSILON;
  Ia *= EPSILON;
  Ra *= EPSILON;
    
  al = -ulEF*an/unEF+(Ra*RuEF+Ia*IuEF)/unEF; // Eq. (3.33) ref [1]
    
  lhs1=(delta*al+2.*sigma*an)/2.;
  lhs2=sin_theta*(SPIN*cos_theta*Ra-r_BL*Ia);

    //std::cout << unEF*al+ulEF*an-RuEF*Ra-IuEF*Ia << std::endl;
  fcov[0] = (lhs1-SPIN*lhs2)/sigma;
  fcov[1]=al/2.-sigma*an/delta;
    
    
  /* CHECK SIGN */
  fcov[2]=(r_BL*Ra+SPIN*cos_theta*Ia);
  fcov[3]=-(SPIN*SQR(sin_theta)*lhs1-r2pa2*lhs2)/sigma;

  // Contravariant SF components
  f[0] = (gtphi*fcov[3]-gphiphi*fcov[0])/detg;
  f[1] = delta*fcov[1]/sigma;
  f[2] = fcov[2]/sigma;
  f[3] = (gtphi*fcov[0]-gtt*fcov[3])/detg;
    
/*---- 
 Evolution of the constants of motion wrt observer's time
       !!!![Check units, by default are  in MBH units]
 ----*/
    
  *dEndt= -fcov[0]/u[0];
  
  *dLzdt = fcov[3]/u[0];
    
  *dQdt = -2.*(LZ-SPIN*E)*((*dLzdt)-SPIN*(*dEndt))
            +2.*((r2pa2*E-LZ*SPIN)*r2pa2*(*dEndt)+(SPIN2*LZ-SPIN*r2pa2*E)*(*dLzdt))/delta
            -2.*fcov[1]*ucov[1]*delta/u[0];
    
/*----- Asigning Global SF values -----*/
  PNFcomp[0] = f[0];  // f_t
  PNFcomp[1] = f[1];  // f_r
  PNFcomp[2] = f[2];  // f_theta
  PNFcomp[3] = f[3];  // f_phi
    
    
/*----- This is the end -----*/
    return 0;
    
}
