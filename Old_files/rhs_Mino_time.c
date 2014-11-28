//
//  Mino_time.c
//  X_ROK_V1.6
//
//  Created by Priscilla on 29/06/2013.
// Based in Drasco and Hughes PRD 73, 024027 (2006) & PRD 69, 044015 (2004)
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"

#include "physical_quantities.h"
#include "macros_2.h"

double ellK(double, double);
double ellE(double, double);

int rhs_Mino_time(double *dlmd_t, double *dlmd_psi, double *dlmd_chi)
{
  double q1,q2,q3,q4,q5,q6,delta,z1;
  double cos_psi,sin_psi,cos_chi;
  double beta,omec2,omEn2;
  double psi, chi;
/*----- Auxiliary quantities -----*/
  psi = PSI_o_NEW;
  chi = CHI_o_NEW;
    
  omEn2 = 1.0 - SQR(ENERGY);    // 1 - E^2
  beta = SQR(SPIN)*omEn2;      // a^2 (1 - E^2)
  omec2 = 1.0 - SQR(ECCENTRICITY);// 1 - e^2
     
  cos_psi = cos(psi);
  sin_psi = sin(psi);
  cos_chi = cos(chi);
  
  q1 = 1.0 + ECCENTRICITY*cos_psi;
  q2 = P_p/q1;           // r = p M / ( 1 + e cos(psi) )
  q3 = P_p - P3_p - ECCENTRICITY*(P_p + P3_p*cos_psi);  // p - p3 - e ( p + p3 cos(psi) )
  q4 = P_p - P4_p + ECCENTRICITY*(P_p - P4_p*cos_psi);  // p - p4 + e ( p - p4 cos(psi) )
  
  q5 = SQR(q2) + SQR(SPIN);// r^2 + a^2
  delta = q5 - 2.0*q2;     // Delta = r^2 - 2 M r + a^2
  
  q6 = ENERGY*(SQR(q5)/delta-SQR(SPIN)) - 2.0*SPIN*q2*LZ/delta;  // gamma = E [ ( r^2 + a^2 )^2 / Delta - a^2 ] - 2 M a r Lz / Delta
  
  z1 = Z_minus*SQR(cos_chi); // z = z_minus cos(chi)^2 = cos(theta)^2
  
/*----- dlambda/dt -----*/
  *dlmd_t = 1.0/(q6 + SQR(SPIN)*ENERGY*z1);
  
/*----- dlamda/dpsi -----*/
 *dlmd_psi = omec2/sqrt(omEn2*q3*q4);
  
/*----- dlamda/dchi -----*/
 *dlmd_chi = sqrt(beta*(Z_plus-z1));
  
/*----- This is the end -----*/
  return 0;
}




//========== OLD ROUTINE=======


/*
double in_psi;        // integration limits
double psi_aux,h;         // Used to keep psi between 0 and 2Pi
double c1 = sqrt(1.- SQR(E));
double c2 = 1.- SQR(e);
double dLr0,dLrF,dLrI,dLraux;// dpsi/dpsidlambda
double int_points = 10;    // integration points,
double deltaPsi;       // increment of psi
int interval_points;     // points in the integration interval
int k;
//int loop =0;        // evolution counter


//-- Computing Lamba_r
deltaPsi = 2.*PI/(psi + 0.001);  // psi_max/psi  the 0.001 factor is for aboinding problems at psi=0
interval_points = ceil(int_points/deltaPsi);

in_psi=0;
dLr0 = c2;
dLr0 /= c1*sqrt(((p-p3) - e*(p + p3*cos(in_psi)))*( (p-p4) +e*(p-p4*cos(in_psi))));

in_psi= psi;
dLrF = c2;
dLrF /= c1*sqrt(((p-p3) - e*(p + p3*cos(in_psi)))*( (p-p4) +e*(p-p4*cos(in_psi))));

in_psi = k*psi/int_points;

for(k=1;k<=(int_points -1.);k++)
{
  dLraux = c2;
  dLraux /= c1*sqrt(((p-p3) - e*(p + p3*cos(in_psi)))*( (p-p4) +e*(p-p4*cos(in_psi))));
  
  dLrI += dLraux;
}

*lambda_r = (psi+0.001)/int_points*(  0.5*(dLr0 + dLrF) + dLrI );

// printf("lambda_r=%4.6e dLr0 =%4.6e dLrF=%4.6e dLrI=%4.6e \n", *lambda_r, dLr0 , dLrF, dLrI);

//-- Computing Lamba_theta
double K, F;    // Complete and incomplete Legendre elliptical integrals
double mode_th, kth;   // Arguments of K and F
double Lamba_theta0;
double betha = SQR(a)*(1.-SQR(E));
//double new_phi;

//---Computing Legendre elliptical integrals First kind---
// mode_th = Z_minus/Z_plus;
kth = sqrt(mode_th);

K  = ellK(0.5*PI,kth);

if( (chi>= 0 )&&(chi<= PI/2) )
{
  F = ellK(0.5*PI-chi,kth);
  Lamba_theta0 = k-F;
  Lamba_theta0 /= sqrt(betha*Z_minus);
  
  *lambda_theta = Lamba_theta0;
  
  // printf("chi=%4.6e, PI=%4.6e  \n", chi,PI/2);
}
else if( (chi>= PI/2)&&(chi<=PI) )
{
  F = ellK(0.5*PI-(PI-chi),kth);
  Lamba_theta0 = k-F;
  Lamba_theta0 /= sqrt(betha*Z_minus);
  
  *lambda_theta = 2*K/sqrt(betha*Z_minus) - Lamba_theta0;
  
  //  printf("chi=%4.6e, PI=%4.6e  \n", chi,PI);
}
else if( (chi>= PI)&&(chi<=3*PI/2) )
{
  F = ellK(0.5*PI-(-PI+chi),kth);
  Lamba_theta0 = k-F;
  Lamba_theta0 /= sqrt(betha*Z_minus);
  
  *lambda_theta = 2*K/sqrt(betha*Z_minus) - Lamba_theta0;
  
  // printf("chi=%4.6e, 3*PI/2=%4.6e  \n", chi,3*PI/2);
}
else if( (chi>= 3*PI/2)&&(chi<=2*PI) )
{
  F = ellK(0.5*PI-(2*PI-chi),kth);
  Lamba_theta0 = k-F;
  Lamba_theta0 /= sqrt(betha*Z_minus);
  
  *lambda_theta = 4*K/sqrt(betha*Z_minus) - Lamba_theta0;
  
  
  //printf("chi=%4.6e, 3*PI/2=%4.6e  \n", chi, 2*PI);
}
else
{ printf("chi=%4.6e  out of bounds!!\n", chi); exit(0);}

*/





