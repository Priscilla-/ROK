/*=====================================================================================
  Finding the partial derivatives of the constant of motion in terms of the orbital 
  elements using the implicit function theorem
 ======================================================================================
 
 Last Update : 10.08.13 by PCM
 ====================================================================================*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/types.h>
#include<time.h>

#include "global_quantities_kerr.h"
#include "global_quantities_main.h"
#include "global_quantities_osc.h"
#include "physical_quantities.h"
#include "macros_2.h"
#include "macros_1.h"


int turning_points_derivatives(double *drpd , double *drad,double *dr3d,double *dr4d,double *dzmd, double *dzpd, double *dpd, double *decd, double *dincd, double *x)
{
  double p,ecc,cosinc,sininc,E,Q,r3,r4,tansqi;
  double ope,ome,ra,rp,ddr2di_ra,ddr2di_rp,tiltcoef,det;
  double prod,beta,omE2,omE4; 
  double Delta_apo,Delta_peri;
  double dVde_ra,dVde_rp,dVdp_ra,dVdp_rp;
  double ddr2dE_ra,ddr2dE_rp,ddr2dLz_ra,ddr2dLz_rp;
  double dEdp,dEde,dEdi,dLZdp,dLZde,dLZdi,dQdp,dQde,dQdi;
  double dpdE, dpdLz, dpdQ,decdE, decdLz, decdQ, dincdE, dincdLz, dincdQ;
  double dpddE, dpddLz, dpddQ;
  double dsrdE, dsrdLz, dsrdQ;
  double dradE, dradLz, dradQ, drpdE, drpdLz,drpdQ;
  double dr3dE, dr3dLz, dr3dQ,dr4dE, dr4dLz, dr4dQ;
  double zm,zp, dzmdE, dzmdLz, dzmdQ,dzpdzm, dzpdE, dzpdLz, dzpdQ;
  //double dpdt, dedt, dincdt;
  
  E = ENERGY;        // Energy of the Orbit
  Q = Q_CONSTANT;    //LZ*LZ*tansqi;     // Carter constant
  
  ra = R_apo;//p/ome;             // Root radial motion r_max (r_apocenter)
  rp = R_peri;//p/ope;             // Root radial motion r_min (r_pericenter)
  r3 = R_3;//params[5];       // Root radial motion
  r4 = R_4;//params[6];       // Root radial motion
  
  zm = Z_minus;
  zp = Z_plus;
  
  p = x[3];            // Orbital parameter semilatus-rectum
  ecc = x[4];          // Orbital parametereccentricity
  
  cosinc = cos(x[5]);  // Cosinus of the orbital parameter tilt
  sininc = sqrt(1. -SQR(cosinc));
  tansqi = 1./(cosinc*cosinc)-1.;
 
  ome = 1.-ecc;           // auxiliar variable
  ope = 1.+ecc;           // auxiliar variable
  
  omE2 = 1.-E*E;
  beta = SQR(SPIN)*omE2;
  omE4 = omE2*omE2;
  prod = r3*r4;
  
  /*---Delta = r^2-2*r+a^2 evaluated @ pericenter and @ apocenter---*/
  Delta_apo  = ra*ra - 2.*ra + SQR(SPIN);
  Delta_peri = rp*rp - 2.*rp + SQR(SPIN);   
  
  //==== Computing d(E,Lz,Q)/ d(e,p,inc) using the implicit function theorem ======//
  /*---- Derivative of the radial potential V(r)= (1-E^2)(r_1-r)(r-r_2)(r-r3)(r-r4) wrt (e,p) @ pericenter and apocenter----*/
  dVde_ra = (1.-E*E)*(ra-rp)*(ra-r3)*(ra-r4)*p/(ome*ome);
  dVde_rp = (1.-E*E)*(ra-rp)*(rp-r3)*(rp-r4)*p/(ope*ope);
  
  dVdp_ra = dVde_ra*ome/p ;
  dVdp_rp = dVde_rp*ope/p;
  
  /*---- Derivative of the radial geodesic equation (drdlambda)^2 = [E(r^2+a^2)-aL_z]^2-Delta[mur^2+(L_z-aE)^2 +Q]= V(r) wrt (E, Lz,inc) @ pericenter and apocenter----*/
  ddr2dE_ra = 2.*(E*(ra*ra+SQR(SPIN))- SPIN*LZ)*(ra*ra+SQR(SPIN))+2.*SPIN*Delta_apo*( LZ-SPIN*E);
  ddr2dE_rp = 2.*(E*(rp*rp+SQR(SPIN))- SPIN*LZ)*(rp*rp+SQR(SPIN))+2.*SPIN*Delta_peri*(LZ-SPIN*E);
  
  ddr2dLz_ra = -2.*SPIN*(E*(ra*ra+SQR(SPIN))-SPIN*LZ)-2*Delta_apo*( LZ-SPIN*E+LZ*tansqi);
  ddr2dLz_rp = -2.*SPIN*(E*(rp*rp+SQR(SPIN))-SPIN*LZ)-2*Delta_peri*(LZ-SPIN*E+LZ*tansqi);
  
  tiltcoef = LZ*LZ*sqrt(tansqi)/(SQR(cosinc));
  ddr2di_ra = 2.*Delta_apo*tiltcoef;                                    
  ddr2di_rp = 2.*Delta_peri*tiltcoef;                    
  
  //--- Apliying Cramer's rule ---//
  
  det = ddr2dE_ra*ddr2dLz_rp-ddr2dLz_ra*ddr2dE_rp;
  
  /*--- Partial derivatives of the constants of motion ---*/
  dEdp  = (ddr2dLz_rp*dVdp_ra + ddr2dLz_ra*dVdp_rp)/det;
  dLZdp = -(ddr2dE_rp*dVdp_ra + ddr2dE_ra*dVdp_rp)/det;
  
  dEde  = (ddr2dLz_rp*dVde_ra - ddr2dLz_ra*dVde_rp)/det;
  dLZde = (-ddr2dE_rp*dVde_ra + ddr2dE_ra*dVde_rp)/det;
  
  dEdi  = (ddr2dLz_rp*ddr2di_ra - ddr2dLz_ra*ddr2di_rp)/det;
  dLZdi = (-ddr2dE_rp*ddr2di_ra + ddr2dE_ra*ddr2di_rp)/det;
  
  /*---- Derivative of the Carter constant Q = L_z^2tan(inc)^2----*/
  dQdp=2.*Q*dLZdp/LZ;                    
  dQde=2.*Q*dLZde/LZ;
  dQdi=2.*Q*dLZdi/LZ+2.*tiltcoef;
  
  //==== Computing  d(e,p,inc)/d(E,Lz,Q) using dE = dE/dp*dp + dE/de*de + dE/dinc*dinc, etc and applying Cramer's rule ======//
  det = dEdp*(dLZde*dQdi-dLZdi*dQde)-dEde*(dLZdp*dQdi-dLZdi*dQdp)+dEdi*(dLZdp*dQde-dLZde*dQdp);
  
  dpdE  = (dLZde*dQdi-dLZdi*dQde)/det;
  dpdLz = -(dEde*dQdi-dEdi*dQde)/det;
  dpdQ  = (dEde*dLZdi-dEdi*dLZde)/det;
  
  decdE  = -(dLZdp*dQdi-dLZdi*dQdp)/det;
  decdLz = (dEdp*dQdi-dEdi*dQdp)/det;
  decdQ  = -(dEdp*dLZdi-dEdi*dLZdp)/det;
  
  dincdE  = (dLZdp*dQde-dLZde*dQdp)/det;
  dincdLz = -(dEdp*dQde-dEde*dQdp)/det;
  dincdQ  = (dEdp*dLZde-dEde*dLZdp)/det;
  
  printf("%4.6e %4.6e %4.6e %4.6e %4.6e", dincdE,dLZdp,dQde,dLZde,dQdp );
  
  dpddE  = prod*(-2./p + 2.*E*dEdp/omE2+dQdp/Q);
  dpddLz = prod*(-2.*ecc/(ope*ome) + 2.*E*dEde/omE2+dQde/Q);
  dpddQ  = prod*( 2.*E*dEdi/omE2 + dQdi/Q);
  
  
  decd[0] = decdE;  decd[1] = decdLz;  decd[2] = decdQ;
  dpd[0]  = dpdE;   dpd[1]  = dpdLz;   dpd[2]  = dpdQ;
  dincd[0]= dincdE; dincd[1]= dincdLz; dincd[2]= dincdQ;

  
  dsrdE  = 4.*E*dEdp/omE4 - 2./(ope*ome);
  dsrdLz = 4.*E*dEde/omE4 - 4.*p*ecc/(SQR(ope*ome));
  dsrdQ  = 4.*E*dEdi/omE4;
  
/*---- radial turning points derivatives----*/                               
  drpdE = dpdE/ope - p*decdE/SQR(ope);
  dradE = dpdE/ome + p*decdE/SQR(ome);
  
  drpdLz = dpdLz/ope - p*decdLz/SQR(ope);
  dradLz = dpdLz/ome + p*decdLz/SQR(ome);
                                 
  drpdQ = dpdQ/ope - p*decdQ/SQR(ope);
  dradQ = dpdQ/ome + p*decdQ/SQR(ome);

  dr3dE = (dpddE - r3*dsrdE)/(r4-r3);
  dr4dE = (dpddE - r4*dsrdE)/(r3-r4);
  
  dr3dLz = (dpddLz - r3*dsrdLz)/(r4-r3);
  dr4dLz = (dpddLz - r4*dsrdLz)/(r3-r4);
  
  dr3dQ = (dpddQ - r3*dsrdQ)/(r4-r3);
  dr4dQ = (dpddQ - r4*dsrdQ)/(r3-r4);
                                 
//==== Computing  d(z_, z+)/d(E,Lz,Q) using implicit function theorem  ======//                                 
/*---- polar turning points derivatives----*/
                                 
  dzmdE = -SQR(SPIN)*E*zm*(1.-zm);
  dzmdE/= beta*(zp-zm);
                                 
  dzmdLz = 2.*LZ*zm;
  dzmdLz/= beta*(zp-zm);

  dzmdQ = -(1.-zm);
  dzmdQ = beta*(zp-zm);
   
  dzpdzm = -Q/(beta*zm*zm);
 
  dzpdE  = dzpdzm*dzmdE;
  dzpdLz = dzpdzm*dzmdLz;
  dzpdQ  = dzpdzm*dzmdQ;

  drad[0] = dradE;  drad[1] = dradLz;  drad[2] = dradQ;
  drpd[0] = drpdE;  drpd[1] = drpdLz;  drpd[2] = drpdQ;
  dr3d[0] = dr3dE;  dr3d[1] = dr3dLz;  dr3d[2] = dr3dQ;
  dr4d[0] = dr4dE;  dr4d[1] = dr4dLz;  dr4d[2] = dr4dQ;
  dzmd[0] = dzmdE;  dzmd[1] = dzmdLz;  dzmd[2] = dzmdQ;
  dzpd[0] = dzpdE;  dzpd[1] = dzpdLz;  dzpd[2] = dzpdQ;
    
      
/*----- this is the end -----*/
return 0;
}
