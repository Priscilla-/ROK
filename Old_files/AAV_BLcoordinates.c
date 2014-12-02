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
 Last Update : 27/04/13
 =============================================================================*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/types.h>
#include<time.h>

#include "global_quantities_kerr.h"
#include "global_quantities_resonances.h"
#include "global_quantities_main.h"
#include "macros_2.h"

void JacobianEllipticFunctions(double, double, double *, double * , double * );

int AAV_BLcoordinates(double W_r, double W_th, double *r_p, double *costh_p)
{
   // FILE *data;
    
    int sign_r, sign_th;
    double varphi_r, varphi_th;
    double sn_r, sn_th;                     // Jacobi's elliptic functions
    double cn_r, cn_th;                     // Jacobi's elliptic functions
    double dn_r, dn_th;                     // Jacobi's elliptic functions
     
/*---Arguments of the Jacobi elliptic functions----*/
    if( (W_r>= 0)&&(W_r<=PI) )
    {
        varphi_r = W_r*Kcom_r/PI;
        sign_r = 1;
    }
    else if( (W_r>= PI)&&(W_r<=2*PI) )
    {
        varphi_r = (2*PI-W_r)*Kcom_r/PI;
        sign_r = -1;
    }
    
    if( (W_th>= 0)&&(W_th<=PI/2) )
    {
        varphi_th = W_th*2.*Kcom_th/PI;
        sign_th = 1;
    }
    else if( (W_th>= PI/2)&&(W_th<=3*PI/2) )
    {
        varphi_th = (PI-W_th)*2.*Kcom_th/PI;
        sign_th =-1;
    }
    else if( (W_th>= 3*PI/2)&&(W_th<=PI) )
    {
        varphi_th = (W_th-2*PI)*2*Kcom_th/PI;
        sign_th = 1;
    }
    
/*---Computing Jacobi elliptic functions----*/
    JacobianEllipticFunctions(varphi_r, SQR(KR),  &sn_r,  &cn_r ,  &dn_r);
    JacobianEllipticFunctions(varphi_th,SQR(KTH), &sn_th, &cn_th , &dn_th);
    
/*----- Boyer-Lindquist radial and polar coordinates ---*/
    *r_p = SQR(sn_r)*(R_apo - R_peri)*R_3 - (R_apo - R_3)*R_peri;
    *r_p  /= SQR(sn_r)*(R_apo - R_peri)-(R_apo - R_3);
    
    *costh_p = sqrt(Z_minus)*sn_th;
        
/*-----This is the end -----*/
  return 0;

}

 