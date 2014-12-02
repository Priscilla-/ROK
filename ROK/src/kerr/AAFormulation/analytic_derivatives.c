
/*==============================================================================
 This routine computes analytical derivatives [see analytic_r_costh.c]
 of bound timelike geodesics in  Kerr in function of the "Mino" time λ and the 
 constant of motion: E, Lz,Q
 -------------------------------------------------------------------------------
  Created by Priscilla on 18/04/2013
 -------------------------------------------------------------------------------
 Last Update : 14/09/13
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


int analytic_derivatives(double *d2rdlE,double *d2rdlLz,double *d2rdlQ,double *d2costhdlE,double *d2costhdlQ,double *d2costhdlLz,
                         double *drpd, double *drad, double *dr3d, double *dr4d, double *dzmd, double *dzpd,double *dOmRd, double *dOmTHd,
                         double *drdra, double *drdrp, double *drdr3, double *drdr4,double *dcosthdzp,double *dcosthdzm,
                         double *drdw,double *dcosthdw,
                         double Wr, double Wth)
{
  double drpdE, dradE, drpdLz, dradLz, drpdQ, dradQ; // "
  double dr3dE, dr4dE, dr3dLz, dr4dLz, dr3dQ, dr4dQ; // "
  double dzmdE, dzpdE, dzmdLz, dzpdLz, dzmdQ, dzpdQ; // "
  //double drdra, drdrp, drdr3, drdr4,dcosthdzp,dcosthdzm; // Radial and polar derivatives wrt  turning points
  double dkrdra, dkrdrp, dkrdr3, dkrdr4, dkthdzm,dkthdzp;// Derivatives of the argument of the eliptical function
  double omE2,beta, SomE2;                  // Auxiliar varialbles
  double dKr,dKth;                          // Derivatives complete elliptical functions
  double dvpr, dvpz;                        // Derivatives argument Jacobi elliptical function
  double sn_r,sn_th,cn_r,cn_th,dn_r,dn_th;  // Jacobi's elliptic functions
  double A, B, dAdra,dAdrp,dAdr3, dAdr4,dBdra,dBdrp,dBdr3, dBdr4;// Aux variables
  double g1,N,D,dN,dNdra, dNdrp, dNdr3, dNdr4; // Aux variables
  double varphi_r,varphi_th;                // Argument Jacobi's elliptic functions
  double dDdE, dDdLz, dDdQ;                 // Aux variables
  double dOmegaRdE,dOmegaRdLz,dOmegaRdQ;    // Derivatives of the radial fundamental frequencies wrt constants of the motion
  double dOmegaTHdE, dOmegaTHdLz, dOmegaTHdQ;//  "                 polar  "
  double denom,N_wr,N_wth;                  // Aux variables
  int sign_r, sign_th;                      // Aux variables

  omE2 = (1.0-SQR(ENERGY));
  beta = SQR(SPIN)*omE2;
  SomE2 = SQR(omE2);
  
  drpdE = drpd[0]; dradE = drad[0]; drpdLz = drpd[1]; dradLz = drad[1]; drpdQ = drpd[2]; dradQ = drad[2];
  dr3dE = dr3d[0]; dr4dE = dr4d[0]; dr3dLz = dr3d[1]; dr4dLz = dr4d[1]; dr3dQ = dr3d[2]; dr4dQ = dr4d[2];
  dzmdE = dzmd[0]; dzpdE = dzpd[0]; dzmdLz = dzmd[1]; dzpdLz = dzpd[1]; dzmdQ = dzmd[2]; dzpdQ = dzpd[2];

//---Computing derivatives of the arguments of the Elliptic functions (here kr--> 1-kr^2 and kth -->1-kth^2)
  dkrdra = -(R_3-R_4)*(R_3 - R_peri);
  dkrdra /=2.*KR*(R_peri-R_4)*SQR(R_apo - R_3);
  
  dkrdrp = (R_3-R_4)*(R_4-R_apo);
  dkrdrp /=2.*KR*(R_apo - R_3)*SQR(R_peri -R_4);
  
  dkrdr3 =  (R_apo - R_peri)*(R_apo-R_4);
  dkrdr3 /= 2.0*KR*(R_peri-R_4)*SQR(R_apo - R_3);
  
  dkrdr4 = (R_apo - R_peri)*(R_3-R_peri);
  dkrdr4 /= 2.0*KR*(R_apo - R_3)*SQR(R_peri-R_4);
   
  dkthdzm = 1./(2.*KTH*Z_plus);
 
  dkthdzp =-Z_minus/(2.0*KTH*SQR(Z_plus));
  
//---- Derivatives of the Complete elliptical integrals of first kind K(k_r), K(k_th)
  dKr  = Ecom_r/(KR*(  1.0-SQR(KR))) - Kcom_r/KR;
  dKth = Ecom_th/(KTH*(1.0-SQR(KTH)))- Kcom_th/KTH;

//---- Derivatives of the fundamental frequencies  wrt constants od motion
  dOmegaRdE = PI*PI*(-2.*ENERGY*(R_apo-R_3)*(R_peri-R_4) + omE2*( (R_peri-R_4)*(dradE-dr3dE) + (R_apo-R_3)*(drpdE-dr4dE)));
  dOmegaRdE /= OMEGA_r*8.*SQR(Kcom_r);
  dOmegaRdE += -OMEGA_r*dKr*(dkrdra*dradE + dkrdrp*drpdE + dkrdr3*dr3dE + dkrdr4*dr4dE)/Kcom_r;

  dOmegaRdLz = PI*PI*(-2.*ENERGY*(R_apo-R_3)*(R_peri-R_4) + omE2*( (R_peri-R_4)*(dradLz-dr3dLz) + (R_apo-R_3)*(drpdLz-dr4dLz)));
  dOmegaRdLz /= OMEGA_r*8.*SQR(Kcom_r);
  dOmegaRdLz += -OMEGA_r*dKr*(dkrdra*dradLz + dkrdrp*drpdLz + dkrdr3*dr3dLz + dkrdr4*dr4dLz)/Kcom_r;

  dOmegaRdQ = PI*PI*(-2.*ENERGY*(R_apo-R_3)*(R_peri-R_4) + omE2*( (R_peri-R_4)*(dradQ-dr3dQ) + (R_apo-R_3)*(drpdQ-dr4dQ)));
  dOmegaRdQ /= OMEGA_r*8.*SQR(Kcom_r);
  dOmegaRdQ += -OMEGA_r*dKr*(dkrdra*dradQ + dkrdrp*drpdQ + dkrdr3*dr3dQ + dkrdr4*dr4dQ)/Kcom_r;

  dOmegaTHdE =  PI*PI*(dzpdE*beta - Z_plus*(2.0*SPIN2*ENERGY));
  dOmegaTHdE /= OMEGA_th*8.*SQR(Kcom_th);
  dOmegaTHdE += -OMEGA_th*dKth*(dkthdzp*dzpdE + dkthdzm*dzmdE)/Kcom_th;

  dOmegaTHdQ =  PI*PI*dzpdQ*beta;
  dOmegaTHdQ /= OMEGA_th*8.*SQR(Kcom_th);
  dOmegaTHdQ += -OMEGA_th*dKth*(dkthdzp*dzpdQ + dkthdzm*dzmdQ)/Kcom_th;

  dOmegaTHdLz =  PI*PI*dzpdLz*beta;
  dOmegaTHdLz /= OMEGA_th*8.*SQR(Kcom_th);
  dOmegaTHdLz += -OMEGA_th*dKth*(dkthdzp*dzpdLz + dkthdzm*dzmdLz)/Kcom_th;

  dOmRd[0] = dOmegaRdE; dOmRd[1] = dOmegaRdLz; dOmRd[2] = dOmegaRdQ;
  dOmTHd[0]= dOmegaTHdE; dOmTHd[1] = dOmegaTHdLz; dOmTHd[2] = dOmegaTHdQ;
    
/*--- Arguments of the Jacobi elliptic functions sn(varphi,k), cn(varphi,k), dn(varphi,k) ----*/
  
/*----- Adjusting angle variables to the range [0, 2*PI]-----*/
 

 // Wr  = fabs(Wr);
 // Wth = fabs(Wth);
  
  if( Wr > 2.*PI)
  {
    N_wr = floor(Wr/(2.*PI));
    Wr -= N_wr*(2.*PI);
  }
  
  if( Wth > 2.*PI)
  {
    N_wth = floor(Wth/(2.*PI));
    Wth -= N_wth*(2.*PI);
  }

//---Radial argument
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
  // else
  //{    printf("Wr =%4,6e out of bounds!!\n", Wr); exit(0);}

/*---Radial Jacobi elliptic function---*/
  JacobianEllipticFunctions(varphi_r, 1.-KR*KR,  &sn_r,  &cn_r ,  &dn_r);

//---Polar argument
  if( (Wth>= 0)&&(Wth<=0.5*PI) )
  {
    varphi_th = Wth*2.0*Kcom_th/PI;
      sign_th = 1;
  }
  else if( (Wth>= 0.5*PI)&&(Wth<=1.5*PI) )
  {
    varphi_th = (PI-Wth)*2.*Kcom_th/PI;
    sign_th =-1;
  }
  else if( (Wth>= 1.5*PI)&&(Wth<=2.*PI) )
  {
    varphi_th = (Wth-2.*PI)*2.*Kcom_th/PI;
    sign_th = 1;
  }
  // else
  //{
  //printf("Wth =%4,6e out of bounds!!\n",Wth);
  //exit(0);
  //}

/*---Polar Jacobi elliptic functions----*/
  JacobianEllipticFunctions(varphi_th,1.-KTH*KTH, &sn_th, &cn_th , &dn_th);
    
//--- Derivative of r wrt Wr ----
    denom = SQR(sn_r)*(R_apo - R_peri)-(R_apo - R_3);
    
    *drdw = 2.*(sign_r)*sn_r*cn_r*dn_r*(R_apo - R_peri)*(R_apo - R_3)*(R_peri - R_3)*Kcom_r;
    *drdw /= SQR(denom)*PI;
    
//--- Derivative of cos(theta) wrt Wth ---
    
    *dcosthdw = 2.*(sign_th)*sqrt(Z_minus)*cn_th*dn_th*Kcom_th;
    *dcosthdw /= PI;


/*---- Derivatives of the arguments of the Jacobi's functions sn(varphi,k), cn(varphi,k), dn(varphi,k) wrt constants of motion---*/
    
//---Radial argument
  if( (Wr>= 0.)&&(Wr<=PI) )
   dvpr = Wr*dKr/PI;
  else if( (Wr>= PI)&&(Wr<=2.*PI) )
   dvpr = (2.0*PI-Wr)*dKr/PI;
  else
  {
    printf("dvp_r out of bounds!!\n");
    exit(0);
  }

//---Polar argument
  if( (Wth>= 0.)&&(Wth<=0.5*PI) )
   dvpz = Wth*2.0*dKth/PI;
  else if( (Wth>= 0.5*PI)&&(Wth<=1.5*PI) )
   dvpz = (PI-Wth)*2.0*dKth/PI;
  else if( (Wth>= 1.5*PI)&&(Wth<=2.*PI) )
   dvpz = (Wth-2.*PI)*2.*dKth/PI;
  else
  {
    printf("dvp_th out of bounds!!\n");exit(0);
  }


/*---- Partial derivatives of r(w_b)  [r = A/B ] wrt the turning points r3, r4----*/
  A = R_3*(R_apo - R_peri)*sn_r*sn_r - R_peri*(R_apo - R_3);
  B = (R_apo - R_peri)*sn_r*sn_r - (R_apo - R_3);

  dAdra =  R_3*sn_r*sn_r + R_3*(R_apo - R_peri)*2.0*sn_r*cn_r*dn_r*dvpr*dkrdra -R_peri;
  dAdrp = -R_3*sn_r*sn_r + R_3*(R_apo - R_peri)*2.0*sn_r*cn_r*dn_r*dvpr*dkrdrp -R_apo;
  dAdr3 = (R_apo - R_peri)*(sn_r*sn_r + R_3*2.0*sn_r*cn_r*dn_r*dvpr*dkrdr3) + R_peri;
  dAdr4 =  R_3*(R_apo - R_peri)*2.*sn_r*cn_r*dn_r*dvpr*dkrdr4;

  dBdra = (R_apo - R_peri)*2.0*sn_r*cn_r*dn_r*dvpr*dkrdra + sn_r*sn_r -1.;
  dBdrp = (R_apo - R_peri)*2.0*sn_r*cn_r*dn_r*dvpr*dkrdrp - sn_r*sn_r;
  dBdr3 = (R_apo - R_peri)*2.0*sn_r*cn_r*dn_r*dvpr*dkrdr3 + 1.0;
  dBdr4 = (R_apo - R_peri)*2.0*sn_r*cn_r*dn_r*dvpr*dkrdr4;

  *drdra = B*dAdra-A*dBdra;
  *drdra /= B*B; 

  *drdrp = B*dAdrp-A*dBdrp;
  *drdrp /=B*B;

  *drdr3 = B*dAdr3-A*dBdr3;
  *drdr3 /=B*B;

  *drdr4 = B*dAdr4-A*dBdr4;
  *drdr4 /=B*B;

/*---- Partial derivatives of cos(theta)(Wth) wrt turning points----*/
  // dcosthdzp = sqrt(Z_minus)*cn_th*dn_th*dvpz*dkthdzp;
  //dcosthdzm = sqrt(Z_minus)*cn_th*dn_th*dvpz*dkthdzm + 0.5*sn_th/sqrt(Z_minus);

  *dcosthdzp = sqrt(Z_minus)*cn_th*dn_th*dvpz*dkthdzp;
  *dcosthdzm = sqrt(Z_minus)*cn_th*dn_th*dvpz*dkthdzm + 0.5*sn_th/sqrt(Z_minus);

 // printf("%4.6e\n",dvpz ); exit(0);
  
/*---- Computing the second derivatives of the geodesic equations [dr/dλ = N*g1*D/B^2, dcos(theta)[Wth]/dλ] wrt the constant of motion---*/

  g1 = (R_peri-R_3)*(R_apo-R_3)*(R_apo-R_peri);
  N  = sn_r*cn_r*dn_r;
  dN = SQR(cn_r*dn_r) - SQR(sn_r)*( SQR(dn_r) + SQR(KR*cn_r) );
  D  = Kcom_r*OMEGA_r/PI;

  dNdra = dN*dvpr*dkrdra*g1 - N*(R_peri-R_3)*(2.*R_apo-R_peri-R_3);
  dNdrp = dN*dvpr*dkrdrp*g1 - N*(R_apo-R_3)*(R_apo-2.*R_peri-R_3);
  dNdr3 = dN*dvpr*dkrdr3*g1 - N*(R_apo-R_peri)*(R_apo+R_peri-2.*R_3);
  dNdr4 = dN*dvpr*dkrdr4*g1;

  dDdE = ( dKr*(dkrdra*dradE + dkrdrp*drpdE + dkrdr3*dr3dE + dkrdr4*dr4dE)*OMEGA_r + Kcom_r*dOmegaRdE )/PI;
  dDdLz= ( dKr*(dkrdra*dradLz + dkrdrp*drpdLz + dkrdr3*dr3dLz + dkrdr4*dr4dLz)*OMEGA_r + Kcom_r*dOmegaRdLz )/PI;
  dDdQ = ( dKr*(dkrdra*dradQ + dkrdrp*drpdQ + dkrdr3*dr3dQ + dkrdr4*dr4dQ)*OMEGA_r + Kcom_r*dOmegaRdQ )/PI;


  *d2rdlE  = 2.0*(sign_r)*(((dNdra*dradE + dNdrp*drpdE + dNdr3*dr3dE + dNdr4*dr4dE)*D + N*g1*dDdE)*B*B
                      -N*g1*D*2.0*B*(dBdra*dradE + dBdrp*drpdE + dBdr3*dr3dE + dBdr4*dr4dE))/(SQR(B*B));

  *d2rdlLz  = 2.0*(sign_r)*(((dNdra*dradLz + dNdrp*drpdLz + dNdr3*dr3dLz + dNdr4*dr4dLz)*D + N*g1*dDdLz)*B*B
                       -N*g1*D*2.0*B*(dBdra*dradLz + dBdrp*drpdLz + dBdr3*dr3dLz + dBdr4*dr4dLz))/(SQR(B*B));


  *d2rdlQ  = 2.0*(sign_r)*(((dNdra*dradQ + dNdrp*drpdQ + dNdr3*dr3dQ + dNdr4*dr4dQ)*D + N*g1*dDdQ)*B*B
                      -N*g1*D*2.0*B*(dBdra*dradQ + dBdrp*drpdQ + dBdr3*dr3dQ + dBdr4*dr4dQ))/(SQR(B*B));


  *d2costhdlE = 2.0*(sign_th)*( 0.5*dzmdE*cn_th*dn_th*Kcom_th*OMEGA_th/sqrt(Z_minus)
                          +sqrt(Z_minus)*( (-(sn_th*dn_r*dn_r + SQR(KTH*cn_th)*sn_th)*dvpz*(dkthdzp*dzpdE+dkthdzm*dzmdE)*Kcom_th*OMEGA_th
                                            +cn_th*dn_th*(dKth*(dkthdzm*dzmdE+dkthdzp*dzpdE)*OMEGA_th
                                                          +Kcom_th*dOmegaTHdE))))/PI;

  *d2costhdlLz = 2.0*(sign_th)*( 0.5*dzmdLz*cn_th*dn_th*Kcom_th*OMEGA_th/sqrt(Z_minus)
                           +sqrt(Z_minus)*( (-(sn_th*dn_r*dn_r + SQR(KTH*cn_th)*sn_th)*dvpz*(dkthdzp*dzpdLz+dkthdzm*dzmdLz)*Kcom_th*OMEGA_th
                                             +cn_th*dn_th*(dKth*(dkthdzm*dzmdLz+dkthdzp*dzpdLz)*OMEGA_th
                                                           +Kcom_th*dOmegaTHdLz))))/PI;

  *d2costhdlQ = 2.0*(sign_th)*( 0.5*dzmdQ*cn_th*dn_th*Kcom_th*OMEGA_th/sqrt(Z_minus)
                          +sqrt(Z_minus)*( (-(sn_th*dn_r*dn_r + SQR(KTH*cn_th)*sn_th)*dvpz*(dkthdzp*dzpdQ+dkthdzm*dzmdQ)*Kcom_th*OMEGA_th
                                            +cn_th*dn_th*(dKth*(dkthdzm*dzmdQ+dkthdzp*dzpdQ)*OMEGA_th
                                                          +Kcom_th*dOmegaTHdQ))))/PI;

/*----- This is the end -----*/
  return 0;
}

