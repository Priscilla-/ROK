/*==============================================================================
 
  This routine calculates the turning points of the orbit. This has been made 
  with the help of Maple.
 
  NOTE: see below for more information

-------------------------------------------------------------------------------
  Last Update : 31.12.12
==============================================================================*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "global_quantities_kerr.h"

int turning_points(double p_orbit, double eccentricity, double energy, double angular_momentum_z, double q_carter,double z_minus, double *r_3, double *r_4, double *p_3, double *p_4, double *z_plus)

{
    double paux_1,paux_2,paux_3,paux_4;             // Auxiliary quantities
    double paux_5,paux_6,paux_7,paux_8;             // Auxiliary quantities
    double raux_1,raux_2,raux_3;                    // Auxiliary quantities
    
    
    /*----- Calculation of the turning points R_3 and R_4 -----*/
    /**----- Auxiliary quantities -----**/
    paux_1 = 1.0 - eccentricity*eccentricity;                    // 1 - e^2
    paux_2 = 1.0 - energy*energy;                                // 1 - E^2
    paux_3 = 1.0/paux_1;                                         // 1 / ( 1 - e^2 )
    paux_4 = 1.0/paux_2;                                         // 1 / ( 1 - E^2 )
    paux_5 = paux_4 - p_orbit*paux_3;                            // 1 / ( 1 - E^2 ) - p / ( 1 - e^2 )
    paux_6 = p_orbit;                                            // p M^2
    paux_7 = paux_6*paux_6;                                      // p^2 M^4
    paux_8 = SPIN2*paux_1*q_carter*paux_4/paux_7;                // a^2 ( 1 - e^2 ) Q / [ p^2 M^4 ( 1 - E^2 ) ]
    
    /**----- Calculation of R_3 -----**/
    *r_3 = paux_5 + sqrt( paux_5*paux_5 - paux_8 );              // r_3 = M ( 1/(1 - E^2) - p/(1 - e^2) + sqrt[ ( 1/(1 - E^2) - p/(1 - e^2) )^2 - a^2 (1 - E^2) C/(p^2 M^4 (1 - E^2))  ] )
    
    /**----- Calculation of R_4 -----**/
    *r_4 = paux_5 - sqrt( paux_5*paux_5 - paux_8 );              // r_4 = M ( 1/(1 - E^2) - p/(1 - e^2) - sqrt[ ( 1/(1 - E^2) - p/(1 - e^2) )^2 - a^2 (1 - E^2) C/(p^2 M^4 (1 - E^2))  ] )
    
    /**----- Calculation of P_3 -----**/
    *p_3 = (paux_5 + sqrt( paux_5*paux_5 - paux_8 ))*(1.0 - eccentricity); // r_3 = p_3 M/(1 - e) => p_3 = r_3 ( 1 - e ) / M
    
    /**----- Calculation of P_4 -----**/
    *p_4 = (paux_5 - sqrt( paux_5*paux_5 - paux_8 ))*(1.0 + eccentricity); // r_4 = p_4 M/(1 + e) => p_4 = r_4 ( 1 + e ) / M
    
    
    /*----- Calculations of the extrema point Z_plus (Z_minus is already computed in 'pez_to_ELCQ.c') -----*/
    if ( FLAG_Schwarzschild == 'y' )
    {
        raux_2 = 0.0; raux_1 = angular_momentum_z*angular_momentum_z;
        raux_3 = 0.0;    // There is no need to do this, but just to have the variables under control in the compilation
        *z_plus = 0.0;
    }
    else
    {
        /**----- Auxiliary quantities -----**/
        raux_1 = SPIN2*paux_2;                                                      // a^2 ( 1 - E^2 )
        raux_2 = raux_1 + angular_momentum_z*angular_momentum_z + q_carter;         // a^2 ( 1 - E^2 ) + L^2 + Q
        raux_3 = raux_2/raux_1;                                                     // [ a^2 ( 1 - E^2 ) + L^2 + Q ] / [ a^2 (1 - E^2) ]
        
        /**----- Calculation of Z_plus -----**/
        *z_plus = raux_3 - z_minus;
    }
    
/*----- Routine ends Here -----*/
    return 0;
    
}

/*****************************************************************************************/
/*                                                                                       */
/*   The conventions of the ROK Code for the Carter Constant are:                        */
/*                                                                                       */
/*      4 .2       2   2           2                         2   2                       */
/*   rho  r  = [ (r + a )E - a Lz ]  - Delta [ Q + (a E - Lz) + r  ]                     */
/*                                                                                       */
/*   Q = C + ( a E - Lz )^2                                                              */
/*                                                                                       */
/*****************************************************************************************/
/*                                                                                       */
/*   Given the orbital parameters of a geodesic orbit around a Kerr BH,                  */
/*                                                                                       */
/*                                (p,e,theta_inc),                                       */
/*                                                                                       */
/*   where p is the semi-latus rectum, e is the eccentricity, theta_inc is the           */
/*   inclination angle of the orbit, and also the values of the constants of motion      */
/*   (E,Lz,Q/C), this routine computes the values of the turning points of the geodesic  */
/*   equation for the Boyer-Lindquist radial coordinates 'r' and 'theta'.                */
/*                                                                                       */
/*---------------------------------------------------------------------------------------*/
/*                                                                                       */
/*                         SOME TECHNICAL DETAILS                                        */
/*                                                                                       */
/* @ The geodesic equation for the Boyer-Lindquist radial coordinate r is given by:      */
/*                                                                                       */
/*         4   / dr \ 2        2   2         2          2            2                   */
/*      rho   | ---- |   = [E(r + a ) - a Lz] - Delta [r + (Lz - a E)  + C]              */
/*             \ dT /                                                                    */
/*                                                                                       */
/*   where rho and Delta are given by                                                    */
/*                                                                                       */
/*                        2    2    2    2                                               */
/*                     rho  = r  + a  cos (theta)                                        */
/*                                                                                       */
/*                               2            2                                          */
/*                      Delta = r  - 2 M r + a                                           */
/*                                                                                       */
/*   This equation can be rewritten as follows:                                          */
/*                                                                                       */
/*         4   / dr \ 2         2                                                        */
/*      rho   | ---- |   =  (1-E ) (r_apo-r) (r-r_peri) (r-r_3) (r-r_4)                  */
/*             \ dT /                                                                    */
/*                                                                                       */
/*   and where the turning points have been defined to have the following expressions:   */
/*                                                                                       */
/*                                       p M                                             */
/*    Apocenter ->             r_apo = -------                                           */
/*                                      1 - e                                            */
/*                                                                                       */
/*                                       p M                                             */
/*    Pericenter ->           r_peri = -------                                           */
/*                                      1 + e                                            */
/*                                                                                       */
/*                                    p_3 M                                              */
/*    r_3 ->                   r_3 = -------                                             */
/*                                    1 - e                                              */
/*                                                                                       */
/*                                    p_4 M                                              */
/*    r_4 ->                   r_4 = -------                                             */
/*                                    1 + e                                              */
/*                                                                                       */
/*    where p is the semilatus rectum and e is the orbital eccentricity. These           */
/*    parameters are inputs of this geodesic solver, so r_peri and r_apo are given       */
/*    straightforwardly.  Then, we are left with the calculation of p_3 and p_4.  As     */
/*    said above, we assume that we also know the value of the constant of motion        */
/*    (E,Lz,Q/C) [the routine 'pei_to_ELQ.c' computes them].                             */
/*                                                                                       */
/* @ The geodesic equation for the Boyer-Lindquist polar angular coordinate theta is     */
/*   given by:                                                                           */
/*                                                                                       */
/*         4   / d theta \ 2            2          2    2   2            2               */
/*      rho   | --------- |   =  C - cot (theta) Lz  - a cos (theta) (1-E )              */
/*             \    dT   /                                                               */
/*                                                                                       */
/*   It is convenient to perform the following coordinate change:                        */
/*                                                                                       */
/*                                    2                                                  */
/*                             z = cos (theta)                                           */
/*                                                                                       */
/*   In this way, the equation for theta becomes:                                        */
/*                                                                                       */
/*         4   / d theta \ 2      2    2   (z - z_+) (z - z_-)                           */
/*      rho   | --------- |   =  a (1-E ) ---------------------                          */
/*             \    dT   /                        1 - z                                  */
/*                                                                                       */
/*   where the turning points are given by:                                              */
/*                                                                                       */
/*          2    2      2              2    2      2    2     2      2                   */
/*         a (1-E ) + Lz + C + Sqrt[ (a (1-E ) + Lz + C) - 4 a C (1-E ) ]                */
/*  z_- = ---------------------------------------------------------------                */
/*                                 2    2                                                */
/*                                a (1-E )                                               */
/*                                                                                       */
/*                                                                                       */
/*          2    2      2              2    2      2    2     2      2                   */
/*         a (1-E ) + Lz + C - Sqrt[ (a (1-E ) + Lz + C) - 4 a C (1-E ) ]                */
/*  z_+ = ---------------------------------------------------------------                */
/*                                 2    2                                                */
/*                                a (1-E )                                               */
/*                                                                                       */
/*****************************************************************************************/