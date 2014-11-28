/*==============================================================================
 This routine computes the analytical expressions for the Mino time lambda,
 parameterized in different ways: Lamda(r) and Lamba(theta), being r and theta 
 BL coordinates
 
 See Fujita and Hikida Class. Quantum Grav. 26 (2009) 135002
 -------------------------------------------------------------------------------
                                            Created by Priscilla on 24/07/2013
 -------------------------------------------------------------------------------
 Last Update: 28/07/2013                                        Modified by PCM
 =============================================================================*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "physical_quantities.h"
#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "global_quantities_resonances.h"

#include "macros_2.h"

double ellK(double, double);


int Mino_time(double r, double cos_theta, double *LAMBDA)
{
  double F_r, F_th;           // Incomplete elliptic integral of the firs kind
  double y_r, y_th, kr, kth;// Arguments elliptical integrals
  double qE, q13,q24,qL;     // Auxiliar variables
  
 /*--- Mino time parameterized in function or r: Lambda(r) ---*/
  qE = sqrt((1.-SQR(ENERGY))*(R_apo - R_3)*(R_peri - R_4));
  
  y_r = (R_apo - R_3)*(r-R_peri);
  y_r/= (R_apo - R_peri)*(r-R_3);
  y_r = asin(sqrt(y_r));
  
  kr = (R_apo - R_peri)*(R_3 - R_4);
  kr/= (R_apo - R_3)*(R_peri - R_4);
  kr = sqrt(kr);
  
  F_r = ellK(y_r,kr);
  
  LAMBDA[2] = 2.*F_r/qE;
 
/*--- Mino time parameterized in function or theta: Lambda(theta) ---*/
  qL = SPIN*sqrt((1.0 - SQR(ENERGY))*Z_plus);
  
  y_th = cos_theta/sqrt(Z_minus);
  
  kth = sqrt(Z_minus/Z_plus);
  
  F_th = ellK(y_th,kth);
  
  LAMBDA[3] = F_th/qL;
  
  //printf("LAMBDA[r] = %4.6e,  LAMBDA[th]=%4.6e, r = %4.6e, theta =%4.6e\n",LAMBDA[2], LAMBDA[3], r, cos_theta);
  
/*----- This is the end -----*/
  return 0;
}

