/*=======================================================================================
                                                                                       
   Given the constants of motion of a geodesic orbit around a Kerr BH,                 
                                                                                       
                           (E,Lz,C/Q),                                                 
                                                                                       
   where E is Energy, Lz is the Angular Momentum component in the direction of the BH  
   Spin, and C/Q is the Carter constant (note that it is not unique and hence there    
   are several choices for its definitions. Here we use the two most popular ones,     
   they are related by:                                                                
                                                                                       
            Q = C - ( a E - Lz )^2                                                     
                                                                                       
   We ASSUME that both, C and Q, are known.                                            
                                                                                       
   This routine returns the values of the orbital parameters (p,e,theta_inc), where p  
   is the semilatus rectum, e is the eccentricity, and theta_inc is the inclination    
   angle.   It also provides other related information, namely: the turning points of  
   the motion [(r_peri, r_apo) and theta_min], another parametrization of the          
   inclination angle [iota], and the extrema [(r_3,r_4) and theta_plus], which are     
   needed for the ODE solver of the Kerr orbits.  
 
   NOTE: See below for further information
==========================================================================================
                                                                                       
   Last Update : 02.02.12 by PCM                                                              
==========================================================================================*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <gsl/gsl_poly.h>

#include "global_quantities_kerr.h"

int checks_errors(int);

int ELCQ_to_peiz(double energy, double Lz, double c_carter, double q_carter,
                double *p_orbit, double *eccentricity, double *theta_inc, double *inclination,
			    double *r_peri, double *r_apo, double *theta_minus, double *z_minus, 
				double *r_3, double *r_4, double *p_3, double *p_4, double *z_plus)
{
  double signo;                                               // Sign to distinguish Prograde from Retrograde Orbits
  
  double a_0, a_1, a_2, a_3;                                  // Coefficients of the Quartic Polynomial
  double c_alpha, c_beta, c_gamma;                            // Coefficients of the "depressed" Quartic Polynomial 
  double b_0, b_1, b_2;                                       // Coefficients of the Auxiliary Cubic Equation
  double root_c1, root_c2, root_c3;                           // Roots of the Auxiliary Cubic Equation
  double y_1;                                                 // A real Root of the Auxiliary Cubic Equation
  double z_minus_sqrt;
  
  double aux_1, aux_2, aux_3, aux_4, aux_5, aux_6, aux_7;     // Auxiliary Quantities
  double aux_8, aux_9, aux_10, aux_11, aux_12, aux_13;        //     "         "
  double aux_14, aux_15, aux_16, aux_17, aux_18;              //     "         "
  double aux_11a, aux_14a, aux_15a;                           //     "         "
  
  
/*----- Solving for the extrema in the radial motion: (r_peri, r_apo, r_3, r_4) -----*/
/**----- Computation of the coefficients of the quartic equation for the extreme in the radial motion -----**/
  aux_1 = 1.0 - energy*energy;                                // 1 - E^2 
  aux_2 = 1.0/aux_1;                                          // 1 / ( 1 - E^2 )
  aux_3 = SPIN2*aux_1;                                     // a^2 ( 1 - E^2 )
  aux_4 = Lz*Lz;                                              // Lz^2
  
  a_3 = - 2.0*aux_2;                                          // - 2 M / ( 1 - E^2 )
  a_2 = ( aux_3 + aux_4 + q_carter )*aux_2;                   // [ a^2 ( 1 - E^2 ) + Lz^2 + Q ] / ( 1 - E^2 )
  a_1 = - 2.0*c_carter*aux_2;                                 // - 2 M C / ( 1 - E^2 )
  a_0 = SPIN2*q_carter*aux_2;                              // a^2 Q / ( 1 - E^2 )
  
/**----- Computation of the coefficients of the "depressed" quartic associated with the initial quartic -----**/
  aux_5 = a_3*a_3;                                            // ( a_3 )^2
  aux_6 = aux_5*a_3;                                          // ( a_3 )^3
  aux_7 = aux_5*aux_5;                                        // ( a_3 )^4
  
  c_alpha = -0.375*aux_5 + a_2;                               // - (3/8) (a_3)^2     |  0.375 = 3/8
  c_beta = 0.125*aux_6 - 0.5*a_2*a_3 + a_1;                   // (1/8) (a_3)^3 - (1/2) a_2 a_3 + a_1   |  0.125 = 1/8
  c_gamma = -0.01171875*aux_7 + 0.0625*a_2*aux_5 - 0.25*a_1*a_3 + a_0;  // - (3/256) (a_3)^4 + (1/16) (a_3)^2 a_2 - (1/4) a_1 a_3 + a_0   |  0.01171875 = 3/256  0.0625 = 1/16

/**----- Computation of the coefficients of the auxiliary cubic equation -----**/
  aux_8 = c_alpha*c_alpha;                                    // alpha^2
  aux_9 = aux_8*c_alpha;                                      // alpha^3
  aux_10 = c_beta*c_beta;                                     // beta^2

  b_2 = 2.5*c_alpha;                                          // (5/2) alpha
  b_1 = 2.0*aux_8 - c_gamma;                                  // 2 alpha^2 - gamma
  b_0 = 0.5*aux_9 - 0.5*c_alpha*c_gamma - 0.125*aux_10;       // (1/2) alpha^3 - (1/2) alpha gamma - (1/8) beta^2

/**----- Solving the cubic with the GSL -----**/
  gsl_poly_solve_cubic(b_2, b_1, b_0, &root_c1, &root_c2, &root_c3); 
 
/**----- Selecting the first root (the one we are sure is real) -----**/
  y_1 = root_c1;
  
/**----- Computation of the roots of the initial quartic -----**/
  aux_11a = c_alpha + 2.0*y_1;                                // alpha + 2 y_1
  aux_11 = sqrt( aux_11a );                                   // sqrt( alpha + 2 y_1 )
  
  if (aux_11a < 0.0)
    checks_errors(19);
  
  
  aux_12 = 3.0*c_alpha + 2.0*y_1;                             // 3 alpha + 2 y_1
  aux_13 = 2.0*c_beta/aux_11;                                 // 2 beta / sqrt( alpha + 2 y_1 )     
  aux_14a = -aux_12 - aux_13;                                 // - ( 3 alpha + 2 y_1 + 2 beta / sqrt[ alpha + 2 y_1 ] ) 
  aux_15a = -aux_12 + aux_13;                                 // - ( 3 alpha + 2 y_1 - 2 beta / sqrt[ alpha + 2 y_1 ] )                  
						
  if ( (aux_14a < 0.0) || (aux_15a < 0.0) )
    checks_errors(19);							 
								 	 
  aux_14 = sqrt( aux_14a );                                   // sqrt{ - ( 3 alpha + 2 y_1 + 2 beta / sqrt[ alpha + 2 y_1 ] ) }
  aux_15 = sqrt( aux_15a );                                   // sqrt{ - ( 3 alpha + 2 y_1 - 2 beta / sqrt[ alpha + 2 y_1 ] ) }
  
// There are four roots which can be written as follows:
//
//  r_one   = -0.25*a_3 + 0.5*(   aux_11 + aux_14 )           Since aux_11 > 0, aux_14 > 0, and aux_15 > 0, we automatically have the following relations:
//  r_two   = -0.25*a_3 + 0.5*(   aux_11 - aux_14 )               r_one > r_two, r_three > r_four, and r_one > r_four
//  r_three = -0.25*a_3 + 0.5*( - aux_11 + aux_15 )
//  r_four  = -0.25*a_3 + 0.5*( - aux_11 - aux_15 )
  
    
  if ( (2.0*aux_11 + aux_14 - aux_15) >= 0.0 )                // Checking whether r_one is larger than r_three
  {                                                           // In the positive case, r_one > r_three, we have that r_apo = r_one
    *r_apo = -0.25*a_3 + 0.5*( aux_11 + aux_14 );             // - (1/4) a_3 + (1/2) ( sqrt( alpha + 2 y_1 ) + sqrt{ - ( 3 alpha + 2 y_1 + 2 beta / sqrt[ alpha + 2 y_1 ] ) } )  
  
    if ( (2.0*aux_11 - aux_14 - aux_15) >= 0.0 )              //   Now we must check whether r_two is larger than r_three
	{                                                         //   In the positive case, r_peri = r_two, r_3 = r_three, and r_4 = r_four
      *r_peri = -0.25*a_3 + 0.5*( aux_11 - aux_14 );          // - (1/4) a_3 + (1/2) ( sqrt( alpha + 2 y_1 ) - sqrt{ - ( 3 alpha + 2 y_1 + 2 beta / sqrt[ alpha + 2 y_1 ] ) } )  
	  *r_3 = -0.25*a_3 + 0.5*( - aux_11 + aux_15 );           // - (1/4) a_3 + (1/2) ( - sqrt( alpha + 2 y_1 ) + sqrt{ - ( 3 alpha + 2 y_1 - 2 beta / sqrt[ alpha + 2 y_1 ] ) } )  
      *r_4 = -0.25*a_3 + 0.5*( - aux_11 - aux_15 );           //   - (1/4) a_3 + (1/2) ( - sqrt( alpha + 2 y_1 ) - sqrt{ - ( 3 alpha + 2 y_1 - 2 beta / sqrt[ alpha + 2 y_1 ] ) } ) 
    }
    else                                                      //   In the negative case, r_peri = r_three
    {
	  *r_peri = -0.25*a_3 + 0.5*( - aux_11 + aux_15 );        //   - (1/4) a_3 + (1/2) ( - sqrt( alpha + 2 y_1 ) + sqrt{ - ( 3 alpha + 2 y_1 - 2 beta / sqrt[ alpha + 2 y_1 ] ) } )  
	  
	  if ( (2.0*aux_11 - aux_14 + aux_15) >= 0.0 )            //     Now we must check whether r_two is larger than r_four 
	  {                                                       //     In the positive case, r_3 = r_two and r_4 = r_four           
        *r_3 = -0.25*a_3 + 0.5*( aux_11 - aux_14 );           //     - (1/4) a_3 + (1/2) ( sqrt( alpha + 2 y_1 ) - sqrt{ - ( 3 alpha + 2 y_1 + 2 beta / sqrt[ alpha + 2 y_1 ] ) } )  
        *r_4 = -0.25*a_3 + 0.5*( - aux_11 - aux_15 );         //     - (1/4) a_3 + (1/2) ( - sqrt( alpha + 2 y_1 ) - sqrt{ - ( 3 alpha + 2 y_1 - 2 beta / sqrt[ alpha + 2 y_1 ] ) } ) 
      }
	  else 
	  {                                                       //     In the negative case, r_3 = r_four and r_3 = r_two
        *r_3 = -0.25*a_3 + 0.5*( - aux_11 - aux_15 );         //     - (1/4) a_3 + (1/2) ( - sqrt( alpha + 2 y_1 ) - sqrt{ - ( 3 alpha + 2 y_1 - 2 beta / sqrt[ alpha + 2 y_1 ] ) } ) 
        *r_4 = -0.25*a_3 + 0.5*( aux_11 - aux_14 );           //     - (1/4) a_3 + (1/2) ( sqrt( alpha + 2 y_1 ) - sqrt{ - ( 3 alpha + 2 y_1 + 2 beta / sqrt[ alpha + 2 y_1 ] ) } )  
	  }
    }
  }
  else                                                        // and automatically, r_peri = r_one, r_3 = r_two, and r_4 = r_four
  {
    *r_apo = -0.25*a_3 + 0.5*( - aux_11 + aux_15 );
	*r_peri = -0.25*a_3 + 0.5*( aux_11 + aux_14 );
	*r_3 = -0.25*a_3 + 0.5*(   aux_11 - aux_14 );
    *r_4 = -0.25*a_3 + 0.5*( - aux_11 - aux_15 );               
  }


/*----- Solving for the extrema in the polar motion: (z_minus, z_plus) -----*/
  if ( FLAG_Schwarzschild == 'y' )
  {
    *z_minus = 0.0;
	*z_plus = 0.0;
  }
  else
  {
    aux_16 = 0.5*a_2/SPIN2;                                  // [ a^2 ( 1 - E^2 ) + Lz^2 + C ] / [ 2 a^2 ( 1 - E^2 ) ]
    aux_17 = aux_16*aux_16;                                     // { [ a^2 ( 1 - E^2 ) + Lz^2 + C ] / [ 2 a^2 ( 1 - E^2 ) ] }^2
    aux_18 = q_carter*aux_2/SPIN2;                           // C / [ a^2 ( 1 - E^2) ]

    if ( FLAG_equatorial_orbit == 'y' )
      *z_minus = 0.0;
    else
      *z_minus = aux_16 - sqrt( aux_17 - aux_18 ); 
  
    *z_plus = aux_16 + sqrt( aux_17 - aux_18 );
  }


/*----- Computing the orbital parameters -----*/
/**----- Eccentricity -----**/
  if ( FLAG_CircOrb == 'y' )
    *eccentricity = 0.0;
  else
    *eccentricity = ( *r_apo - *r_peri )/( *r_apo + *r_peri );
  
/**----- Semilatus Rectum -----**/
  *p_orbit = 2.0*(*r_apo)*(*r_peri)/( *r_apo + *r_peri );
  
/**----- Calculation of P_3 -----**/  
  *p_3 = (*r_3)*(1.0 - *eccentricity);                           // r_3 = p_3 M/(1 - e) => p_3 = r_3 ( 1 - e ) / M
    
/**----- Calculation of P_4 -----**/  
  *p_4 = (*r_4)*(1.0 + *eccentricity);                           // r_4 = p_4 M/(1 + e) => p_4 = r_4 ( 1 + e ) / M
  
/**----- Inclination Angle: Theta_inc -----**/
  z_minus_sqrt = sqrt(*z_minus);	  
  *theta_minus = acos(z_minus_sqrt);
   
  if ( Lz < 0.0 ) signo = -1.0; else signo = 1.0;
  
  *theta_inc = 0.5*PI - signo*(*theta_minus);
 
/**----- Inclination Angle: iota -----**/
  *inclination = acos( Lz/sqrt( Lz*Lz + q_carter ) );
    
   // printf("%4.6e  %4.6e \n",Lz, cos(*inclination) ); exit(0);
/*----- this is the end -----*/
  return 0;
}

/*****************************************************************************************/
/*                                                                                       */
/*   The conventions of the ROK Code for the Carter Constant are:                      */
/*                                                                                       */
/*      4 .2       2   2           2                         2   2                       */
/*   rho  r  = [ (r + a )E - a Lz ]  - Delta [ Q + (a E - Lz) + r  ]                     */
/*                                                                                       */
/*   Q = C - ( a E - Lz )^2                                                              */
/*                                                                                       */
/*****************************************************************************************/
/*                                                                                       */
/*   Given the constants of motion of a geodesic orbit around a Kerr BH,                 */
/*                                                                                       */
/*                           (E,Lz,C/Q),                                                 */
/*                                                                                       */
/*   where E is Energy, Lz is the Angular Momentum component in the direction of the BH  */
/*   Spin, and C/Q is the Carter constant (note that it is not unique and hence there    */
/*   are several choices for its definitions. Here we use the two most popular ones,     */
/*   they are related by:                                                                */
/*                                                                                       */
/*            Q = C - ( a E - Lz )^2                                                     */
/*                                                                                       */
/*   We ASSUME that both, C and Q, are known.                                            */
/*                                                                                       */
/*   This routine returns the values of the orbital parameters (p,e,theta_inc), where p  */
/*   is the semilatus rectum, e is the eccentricity, and theta_inc is the inclination    */
/*   angle.   It also provides other related information, namely: the turning points of  */
/*   the motion [(r_peri, r_apo) and theta_min], another parametrization of the          */
/*   inclination angle [iota], and the extrema [(r_3,r_4) and theta_plus], which are     */
/*   needed for our ODE solver of the Kerr orbits.  To see how they appear in the        */
/*   geodesic equations, have a look at the beginning of the file:                       */
/*                                                                                       */
/*                        compute_extrema_of_motion.c                                    */
/*                                                                                       */
/*****************************************************************************************/
/*                                                                                       */
/*   Given the constants of motion (E,Lz,C/Q) we can recover the information of the      */
/*   orbital parameters (p,e,theta_inc) and derived quantities by looking at the extrema */
/*   of both the radial and polar motion.  That this, the roots of the (dr/d tau) = 0,   */
/*   and (d theta/d tau) = 0.  In the case of the radial motion, this means to find the  */
/*   four real roots of the following quartic polynomial:                                */
/*                                                                                       */
/*   ( 1 - E^2 ) r^4 - 2 M r^3 + [ a^2 ( 1 - E^2 ) + Lz^2 + C ] r^2 - 2 M Q r + a^2 C    */
/*   = 0 ,                                                                               */
/*                                                                                       */
/*   which can be written as:                                                            */
/*                                                                                       */
/*   r^4 + a_3 r^3 + a_2 r^2 + a_1 r + a_0 = 0 ,                                         */
/*                                                                                       */
/*   where the coefficients are:                                                         */
/*                                                                                       */
/*   a_3 = - 2 M / ( 1 - E^2 ) ,                                                         */
/*                                                                                       */
/*   a_2 = [ a^2 ( 1 - E^2 ) + Lz^2 + C ] / ( 1 - E^2 ) ,                                */
/*                                                                                       */
/*   a_1 = - 2 M Q / ( 1 - E^2 ) ,                                                       */
/*                                                                                       */
/*   a_0 = a^2 C / ( 1 - E^2 ) .                                                         */
/*                                                                                       */
/*   We can solve this quartic equation by steps.  First, we need to consider the        */
/*   following cubic equation:                                                           */
/*                                                                                       */
/*   y^3 + b_2 y^2 + b_1 y + b_0 = 0 ,                                                   */
/*                                                                                       */
/*   where the coefficients are given by                                                 */
/*                                                                                       */
/*   b_2 = (5/2) alpha ,                                                                 */
/*                                                                                       */
/*   b_1 = 2 alpha^2 - gamma ,                                                           */
/*                                                                                       */
/*   b_0 = (1/2) alpha^3 - (1/2) alpha gamma - (1/8) beta^2 ,                            */
/*                                                                                       */
/*   where alpha, beta, and gamma are the coefficients of the "depressed" quartic        */
/*   associated with the initial quartic via the following variable transformation:      */
/*                                                                                       */
/*   r = u - (1/4) a_3 ,                                                                 */
/*                                                                                       */
/*   which yields the "depressed" quartic                                                */
/*                                                                                       */
/*   u^4 + alpha u^2 + beta u + gamma = 0 .                                              */
/*                                                                                       */
/*   Then, the relations between the coefficients of the "depressed" quartic and those   */
/*   of the initial quartic are:                                                         */
/*                                                                                       */
/*   alpha = - (3/8) (a_3)^2 + a_2 ,                                                     */
/*                                                                                       */
/*   beta = (1/8) (a_3)^3 - (1/2) a_2 a_3 + a_1 ,                                        */
/*                                                                                       */
/*   gamma = - (3/256) (a_3)^4 + (1/16) a_2 (a_3)^2 - (1/4) a_1 a_3 + a_0 .              */
/*                                                                                       */
/*   Now, let us consider a real solution the cubic equation above, namely y_1 [Note     */
/*   a cubic equation with real coefficients always has at least a real root].  Then,    */
/*   the four solutions of the initial quartic can be written as follows:                */
/*                                                                                       */
/*                                                                                       */
/*    r = - (1/4) a_3 + (1/2) { s_1 sqrt[ alpha + 2 y_1 ] + s_2 sqrt[ - ( 3 alpha        */
/*      + 2 y_1 + s_1 ( 2 beta ) / sqrt( alpha + 2 y_1 ) ) ] } ,                         */
/*                                                                                       */
/*   where s_1 and s_2 are two independent signs, which leads to the four solutions.     */
/*                                                                                       */
/*                                                                                       */
/*   In the case of the polar motion, things are a bit easier.  We need to find the      */
/*   roots of the following quadratic equation for z = cos^2 (theta):                    */
/*                                                                                       */
/*   z^2 + c_1 z + c_0 = 0 ,                                                             */
/*                                                                                       */
/*   where the two coefficients are given by the following expressions:                  */
/*                                                                                       */
/*   c_1 = - [ a^2 ( 1 - E^2 ) + Lz^2 + C ] / [ a^2 ( 1 - E^2 ) ] ,                      */
/*                                                                                       */
/*   c_0 = C / [ a^2 ( 1 - E^2 ) ] .                                                     */
/*                                                                                       */
/*   The two roots of this equation are:                                                 */
/*                                                                                       */
/*   z_- = - ( c_1 / ( 2 c_0 ) ) - sqrt{ [ c_1 / ( 2 c_0 ) ]^2 - c_0 }                   */
/*                                                                                       */
/*   z_+ = - ( c_1 / ( 2 c_0 ) ) + sqrt{ [ c_1 / ( 2 c_0 ) ]^2 - c_0 }                   */
/*                                                                                       */
/*                                                                                       */
/*   This completes the question of finding the roots of the polynomial equations        */
/*   associated with the extrema of the radial and polar motions.  Now, let us see how   */
/*   to compute the rest of the parameters.                                              */
/*                                                                                       */
/*   - Eccentricity = ( r_apo - r_peri ) / ( r_apo + r_peri )                            */
/*                                                                                       */
/*   - Semilatus Rectum = 2 r_apo r_peri / ( r_apo + r_peri ) / mass                     */
/*                                                                                       */
/*   - Inclination Angles:                                                               */
/*                                                                                       */
/*     Theta_inc = PI/2 - sign(Lz)*acos( sqrt( Z_minus ) )                               */
/*                                                                                       */
/*     Iota = acos( Lz / sqrt( Lz^2 + Q ) )                                              */
/*                                                                                       */
/*****************************************************************************************/