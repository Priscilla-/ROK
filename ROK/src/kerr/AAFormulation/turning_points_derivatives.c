/*=====================================================================================
 In the osculating evolution formulation, we need to know how the turning points evolve
 with time for each new geodesic orbit. This routine finds the partial derivatives of 
 the constant of motion in terms of the orbital elements using the implicit function 
 theorem and the Cramer's rule.
 ======================================================================================
 
 Last Update : 10.08.13 by PCM
 ====================================================================================*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/types.h>
#include<time.h>
#include <math.h> 

#include "global_quantities_kerr.h"
#include "global_quantities_main.h"
#include "global_quantities_osc.h"
#include "physical_quantities.h"
#include "macros_2.h"
//#include "macros_1.h"

int turning_points_derivatives(double *drpd , double *drad,double *dr3d,double *dr4d,double *dzmd, double *dzpd, double *dpd, double *decd, double *dthind, double *x)
{
  double p,ecc,E,Q,r3,r4;
  double ope,ome,ra,rp,ddr2dthinc_rp,ddr2dthinc_ra,det;
  double prod,beta,omE2,omE4;
  double Delta_apo,Delta_peri;
  double dVde_ra,dVde_rp,dVdp_ra,dVdp_rp;
  double ddr2dE_ra,ddr2dE_rp,ddr2dLz_ra,ddr2dLz_rp;
  double dEdp,dEde,dEdthinc,dLZdp,dLZde,dLZdthinc,dQdp,dQde,dQdthinc;
  double dpdE, dpdLz, dpdQ,decdE, decdLz, decdQ;
  double dpddE, dpddLz, dpddQ;
  double dsrdE, dsrdLz, dsrdQ;
  double dradE, dradLz, dradQ, drpdE, drpdLz,drpdQ;
  double dr3dE, dr3dLz, dr3dQ,dr4dE, dr4dLz, dr4dQ;
  double zm,zp, dzmdE, dzmdLz, dzmdQ,dzpdE, dzpdLz, dzpdQ;
  double sinthinc, costhinc,dzmdthinc;
  double dthincdE,dthincdLz, dthincdQ,dQfactor;

    
  E = ENERGY;        // Energy of the Orbit
  Q = Q_CONSTANT;    // Carter constant
  
  ra = R_apo;   // Root radial motion r_max (r_apocenter)
  rp = R_peri;  // Root radial motion r_min (r_pericenter)
  r3 = R_3;     // Root radial motion
  r4 = R_4;     // Root radial motion
  
  zm = Z_minus; // Root polar motion
  zp = Z_plus;  // Root polar motion
  
  p = x[3];      // Orbital parameter semilatus-rectum
  ecc = x[4];    // Orbital parametereccentricity
  
  costhinc = cos(x[5]);               // cos(theta_inc)
  sinthinc = sqrt(1.- SQR(costhinc)); // sin(theta_inc)
  
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
  
  /*---- Derivative of the radial geodesic equation (drdlambda)^2 = [E(r^2+a^2)-aL_z]^2-Delta[mur^2+(L_z-aE)^2 +Q]= V(r) wrt (E, Lz,inc) @ pericenter and apocenter - note Q is taken to be Q=z_{a^2*(1-E^2) + L^2/(1-z_)} ----*/
  
  ddr2dE_ra = 2.*(E*(ra*ra+SQR(SPIN))- SPIN*LZ)*(ra*ra+SQR(SPIN))+ Delta_apo*2.*SPIN*(LZ+SPIN*E*(zm-1.));
  ddr2dE_rp = 2.*(E*(rp*rp+SQR(SPIN))- SPIN*LZ)*(rp*rp+SQR(SPIN))+ Delta_peri*2.*SPIN*(LZ+SPIN*E*(zm-1.));
  
  ddr2dLz_ra = -2.*SPIN*(E*(ra*ra+SQR(SPIN))-SPIN*LZ)-2*Delta_apo*(LZ-SPIN*E+LZ*zm/(1.-zm));// tansqi);
  ddr2dLz_rp = -2.*SPIN*(E*(rp*rp+SQR(SPIN))-SPIN*LZ)-2*Delta_peri*(LZ-SPIN*E+LZ*zm/(1.-zm));// tansqi);
  
  dzmdthinc = 2.*sinthinc*costhinc; // sin(2*thinc)
  dQfactor = (dzmdthinc*(beta+ LZ*LZ/(1.-zm)) + zm*LZ*LZ*dzmdthinc/SQR(1.-zm));//2.*Q*dLZdi/LZ + tiltcoef;

  ddr2dthinc_ra = -Delta_apo*dQfactor;
  ddr2dthinc_rp = -Delta_peri*dQfactor;
    
  //--- Apliying Cramer's rule to find the constant of motion and orbital parameter derivatives---//
  
  det = ddr2dE_ra*ddr2dLz_rp-ddr2dLz_ra*ddr2dE_rp;
  
  /*--- Partial derivatives of the constants of motion ---*/
  dEdp  = (ddr2dLz_rp*dVdp_ra + ddr2dLz_ra*dVdp_rp)/det;
  dLZdp = -(ddr2dE_rp*dVdp_ra + ddr2dE_ra*dVdp_rp)/det;
  
  dEde  = (ddr2dLz_rp*dVde_ra - ddr2dLz_ra*dVde_rp)/det;
  dLZde = (-ddr2dE_rp*dVde_ra + ddr2dE_ra*dVde_rp)/det;
  
  dEdthinc  = (ddr2dLz_rp*ddr2dthinc_ra - ddr2dLz_ra*ddr2dthinc_rp)/det;
  dLZdthinc = (-ddr2dE_rp*ddr2dthinc_ra + ddr2dE_ra*ddr2dthinc_rp)/det;
  
  /*---- Derivatives of the Carter constant Q=z_{a^2*(1-E^2) + L^2/(1-z_)} wrt orbital parameters---*/
  dQdp = zm*(-SPIN*SPIN*2.*E*dEdp + 2.*LZ*dLZdp/(1.-zm));
  dQde = zm*(-SPIN*SPIN*2.*E*dEde + 2.*LZ*dLZde/(1.-zm));
  dQdthinc = dzmdthinc*(beta + LZ*LZ/(1.-zm))+zm*( -SPIN*SPIN*2.*E*dEdthinc + (2.*LZ*dLZdthinc*(1.-zm) + LZ*LZ*dzmdthinc)/SQR(1.-zm));
  
  //==== Computing  d(e,p,inc)/d(E,Lz,Q) using dE = dE/dp*dp + dE/de*de + dE/dinc*dinc, etc and applying Cramer's rule ======//
  det = dEdp*(dLZde*dQdthinc-dLZdthinc*dQde)-dEde*(dLZdp*dQdthinc-dLZdthinc*dQdp)+dEdthinc*(dLZdp*dQde-dLZde*dQdp);
  
  dpdE  = (dLZde*dQdthinc-dLZdthinc*dQde)/det;
  dpdLz = -(dEde*dQdthinc-dEdthinc*dQde)/det;
  dpdQ  = (dEde*dLZdthinc-dEdthinc*dLZde)/det;
  
  decdE  = -(dLZdp*dQdthinc-dLZdthinc*dQdp)/det;
  decdLz = (dEdp*dQdthinc-dEdthinc*dQdp)/det;
  decdQ  = -(dEdp*dLZdthinc-dEdthinc*dLZdp)/det;
  
  dthincdE = (dLZdp*dQde-dLZde*dQdp)/det;
  dthincdLz= -(dEdp*dQde-dEde*dQdp)/det;
  dthincdQ = (dEdp*dLZde-dEde*dLZdp)/det;
  
  decd[0] = decdE;  decd[1] = decdLz;  decd[2] = decdQ;
  dpd[0]  = dpdE;   dpd[1]  = dpdLz;   dpd[2]  = dpdQ;
  dthind[0]= dthincdE; dthind[1]= dthincdLz; dthind[2]= dthincdQ;

  dpddE  = prod*(-2./p + 2.*E*dEdp/omE2+dQdp/Q);
  dpddLz = prod*(-2.*ecc/(ope*ome) + 2.*E*dEde/omE2+dQde/Q);
  dpddQ  = prod*( 2.*E*dEdthinc/omE2 + dQdthinc/Q);
  
  dsrdE  = 4.*E*dEdp/omE4 - 2./(ope*ome);
  dsrdLz = 4.*E*dEde/omE4 - 4.*p*ecc/(SQR(ope*ome));
  dsrdQ  = 4.*E*dEdthinc/omE4;
  
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
                                 
                                
/*---- polar turning points derivatives----*/

//--- Derivatives of z_ = cos(thmin)^2  
  dzmdE = dzmdthinc*dthincdE;
  dzmdQ = dzmdthinc*dthincdQ;
  dzmdLz = dzmdthinc*dthincdLz;
  
//--- Derivatives of Z+ = Q/(a^2*(1-E^2)*z_)
  dzpdE  = -Q*(-2.*E*zm +(1.-E*E)*dzmdE)/(SQR((1.-E*E)*zm)*SPIN*SPIN);
  dzpdQ  = (beta*zm - Q*beta*dzmdQ)/SQR(beta*zm);
  dzpdLz = -Q*dzmdLz/(beta*zm*zm);

  
  drad[0] = dradE;  drad[1] = dradLz;  drad[2] = dradQ;
  drpd[0] = drpdE;  drpd[1] = drpdLz;  drpd[2] = drpdQ;
  dr3d[0] = dr3dE;  dr3d[1] = dr3dLz;  dr3d[2] = dr3dQ;
  dr4d[0] = dr4dE;  dr4d[1] = dr4dLz;  dr4d[2] = dr4dQ;
  dzmd[0] = dzmdE;  dzmd[1] = dzmdLz;  dzmd[2] = dzmdQ;
  dzpd[0] = dzpdE;  dzpd[1] = dzpdLz;  dzpd[2] = dzpdQ;
    
      
/*----- this is the end -----*/
return 0;
}
