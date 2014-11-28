/*==============================================================================
 This routine computes analytical solutions of bound timelike geodesics in  Kerr
 in function of the "Mino" time λ for MBH spin smaller than MBH mass.
 
 The expressions computed here are quoted in Fujita et al (2009) [1], where
 analitical expressions for the radial solution r(W_r) and polar
 cos_theta(W_theta) variables are obtained in fuction of the action angle
 variables: W_r = Omega_r*λ and W_th = Omega_theta*λ
 
 -------------------------------------------------------------------------------
 Created by Priscilla on 18/04/2013
 -------------------------------------------------------------------------------
 Last Update : 29/07/13
 =============================================================================*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/types.h>
#include<time.h>

#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"
#include "global_quantities_main.h"
#include "macros_2.h"

void JacobianEllipticFunctions(double, double, double *, double * , double * );

int analytic_r_costh(double *x, double *r_p, double *costh_p)
{
    int sign_r, sign_th;
    double varphi_r, varphi_th;
    double sn_r, sn_th;                     // Jacobi's elliptic functions
    double cn_r, cn_th;                     // Jacobi's elliptic functions
    double dn_r, dn_th;                     // Jacobi's elliptic functions
    double W_r, W_th;
    double N_wr, N_wth;
  
    if(MAG_SF != 0.)
    {
      W_r  = x[6];
      W_th = x[7];
    }
    else
    {
      W_r  = x[1];
      W_th = x[2];
    }
    
   if( W_r > 2.*PI)
   {
    N_wr = floor(W_r/(2.*PI));
    W_r -= N_wr*(2.*PI);
   }
   if(W_th > 2.*PI)
   {
    N_wth = floor(W_th/(2.*PI));
    W_th -= N_wth*(2.*PI);
   }
  

/*---Arguments of the Jacobi elliptic functions----*/
  
    //Radial conditions
    if( (W_r >= 0.0 )&&(W_r <= PI) )
    {
        varphi_r = W_r*Kcom_r/PI;
        sign_r = 1;
    }
    else if( (W_r>= PI)&&(W_r<=2.0*PI) )
    {
        varphi_r = (2.0*PI-W_r)*Kcom_r/PI;
        sign_r = -1;
    }
    else 
    {
      printf("W_r =%4.6e > 2PI =%4.6e \n",W_r, 2.0*PI);
      exit(0);
    }
    
    
    //Polar conditions
    if( (W_th>= 0.0)&&(W_th<=0.5*PI) )
    {
        varphi_th = W_th*2.0*Kcom_th/PI;
        sign_th = 1;
    }
    else if( (W_th>= 0.5*PI)&&(W_th<=1.5*PI) )
    {
        varphi_th = (PI-W_th)*2.0*Kcom_th/PI;
        sign_th =-1;
    }
    else if( (W_th>= 1.5*PI)&&(W_th<= 2.*PI) )
    {
        varphi_th = (W_th-2.0*PI)*2.0*Kcom_th/PI;
        sign_th = 1;
    }
    else
    { 
      printf("W_th =%4.6e > 2PI =%4.6e \n",W_th, 2.0*PI);
      exit(0);
    }
    
/*---Computing Jacobi elliptic functions----*/
    JacobianEllipticFunctions(varphi_r, 1.- KR*KR,  &sn_r,  &cn_r ,  &dn_r);
    JacobianEllipticFunctions(varphi_th,1.-KTH*KTH, &sn_th, &cn_th , &dn_th);
    
/*----- Boyer-Lindquist radial and polar coordinates ---*/
    *r_p  = SQR(sn_r)*(R_apo - R_peri)*R_3 - (R_apo - R_3)*R_peri;
    *r_p /= SQR(sn_r)*(R_apo - R_peri)-(R_apo - R_3);
    
    *costh_p = sqrt(Z_minus)*sn_th;
    
/*----- Mapping onto the Radial and polar phases ---*/
  if(MAG_SF != 0)
  {
    PSI_NEW = acos((x[3]/(*r_p) - 1.)/x[4]);
    CHI_NEW = acos((*costh_p)/sqrt(Z_minus)); 
  }
  else
  {
    PSI_NEW = acos((P_p/(*r_p) - 1.)/ECCENTRICITY); 
    CHI_NEW = acos((*costh_p)/sqrt(Z_minus));
  }
  
  
/*  if(*r_p == R_peri)
    printf("pericenter passage: W_r=%4.6e, W_th=%4.6e, r_p=%4.6e, costh=%4.6e \n", W_r, W_th, *r_p, *costh_p);
  if(*costh_p == sqrt(Z_minus))
    printf("Z_ passage: W_r=%4.6e, W_th=%4.6e, r_p=%4.6e, costh=%4.6e \n");
  
*/
  // printf("sn=%4.6e\n\n",sn_th);
/*-----This is the end -----*/
    return 0;
    
}

