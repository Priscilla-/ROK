/*==============================================================================

   Given the orbital parameters of a geodesic orbit around a Kerr BH,                  
                                                                                       
                           (p,e,theta_inc),                                            
                                                                                       
   where p is the semi-latus rectum, e is the eccentricity, and theta_inc is a measure 
   of the tilt angle of the orbit.  This is different from iota, which is       
   defined by:                                                                         
                                                                                       
          tan^2(iota) = C/Lz^2 .                                                       
                                                                                       
   The definition of theta_inc is:                                                     
                                                                                       
          theta_inc = Pi/2 - sign(Lz) theta_min,                                     
                                                                                       
   where theta_min is the minimum value of the polar angle theta, that is, it is a   
   turning point of the theta motion.                                                  
                                                                                       
   This routine returns the values of the constants of motion (E,Lz,C/Q).              
 -------------------------------------------------------------------------------
 
   Last Update : 13.07.13  by PCM                                                            
===============================================================================*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "global_quantities_kerr.h"

int checks_errors(int);

int pez_to_ELCQ(char FLAG_circular, double p_orbit, double eccentricity, double theta_inc, double spin,
			    double *tilt, 
			    double *r_peri, double *r_apo, double *theta_min, double *z_minus,
			    double *energy, double *Lz, double *c_carter, double *q_carter)
{
  double signo;
  double cos_theta_min;                                     // Auxiliary variables constructed out of Theta_minimum
  double tani, cosi;                                          // Tangent  and cosinus of the other tilt angle (iota)

  double ra2, rp2, rp_minus_mass, spin2;                             // Auxiliary variables constructed in terms of fundamental quantities
  double energy_square, angular_momentum_square;              // idem
  double one_minus_z_m, two_mass_ra, two_mass_rp;             // idem
  double delta_a, delta_p, alpha_a, alpha_p;                  // idem
  double beta_a, beta_p, gamma_a, gamma_p;                    // idem
  double lambda_a, lambda_p, alpha_cross_beta;                // idem 
  double alpha_cross_gamma, lambda_cross_beta;                // idem
  double lambda_cross_gamma, gamma_cross_beta;                // idem
  
  double numer1, numer2, denom1;                              // Other auxiliary quantities
  
  spin2 = spin*spin;
  
/*----- Setting the value of the minimal Polar Angle (Theta_minus) and of z_ (Z_minus) -----*/
  if ( FLAG_equatorial_orbit == 'y' )
  {
    *theta_min = 0.5*PI;
	  cos_theta_min = 0.0;
	 
	  *z_minus = 0.0;
  }
  else
  {
    if ( theta_inc >= 0.0 ) 
	  signo = 1.0;
    else if (theta_inc < 0.0 ) 
	  signo = -1.0;

    *theta_min = 0.5*PI - signo*theta_inc;
    cos_theta_min = cos(*theta_min);
  
    *z_minus = cos_theta_min*cos_theta_min;
  }

/*----- Circular Orbit Case (r = const.) -----*/
  if ( FLAG_circular == 'y' )
  {
/**----- Setting the value of the radius of the circular orbit -----*/
    *r_peri = p_orbit;                                       // For a circular orbit r_peri = r_apo  
    *r_apo = p_orbit;                                        //  "  "    "       "      "   "    "
 
/**----- Constructing Auxiliary Quantities -----**/
    rp2 = (*r_peri)*(*r_peri);                               // r_peri^2
    one_minus_z_m = 1.0 - *z_minus;                          // 1 - z-
    two_mass_rp = 2.0*(*r_peri);                             // 2 M r_peri
    delta_p = rp2 + spin2 - two_mass_rp;                     // r_peri^2 - 2 M r_peri + a^2
	  rp_minus_mass = *r_peri - 1.0;                           // r_peri - M
  
/**----- Constructing the coefficients of the equations for (E,Lz) [after substituting C(E,Lz)] -----**/
    alpha_a = (rp2 + spin2)*(rp2 + spin2*(*z_minus)) + two_mass_rp*spin2*one_minus_z_m;
    alpha_p = 2.0*(*r_peri)*(rp2 + spin2) - spin2*rp_minus_mass*one_minus_z_m;
  
    beta_a = - two_mass_rp*spin;
    beta_p = - spin;
  
    gamma_a = - ( rp2 + spin2*(*z_minus) - two_mass_rp );
    gamma_p = - rp_minus_mass;
  
    lambda_a = - ( rp2 + spin2*(*z_minus) )*delta_p;
    lambda_p = - (*r_peri)*delta_p - rp_minus_mass*( rp2 + spin2*(*z_minus) );
  
/**----- 'Cross Products': Coefficients of the combined equations: A_cross_B = A_a B_p - A_p B_a -----**/  
    alpha_cross_gamma = alpha_a*gamma_p - alpha_p*gamma_a;
    alpha_cross_beta = alpha_a*beta_p - alpha_p*beta_a;
    gamma_cross_beta = gamma_a*beta_p - gamma_p*beta_a;
    lambda_cross_gamma = lambda_a*gamma_p - lambda_p*gamma_a;
    lambda_cross_beta = lambda_a*beta_p - lambda_p*beta_a;
  
/**----- Computing the Energy -----**/
/***----- Auxiliary Quantities: E^2 = ( numer1 +- 2 sqrt(numer2) )/denom1 -----***/
    numer1 = - ( alpha_cross_gamma*lambda_cross_gamma + 2.0*one_minus_z_m*gamma_cross_beta*lambda_cross_beta );
  
    numer2 = one_minus_z_m*gamma_cross_beta*( lambda_cross_beta*( alpha_cross_gamma*lambda_cross_gamma 
	       + one_minus_z_m*gamma_cross_beta*lambda_cross_beta ) - alpha_cross_beta*lambda_cross_gamma*lambda_cross_gamma );
		 
    denom1 = alpha_cross_gamma*alpha_cross_gamma + 4.0*one_minus_z_m*alpha_cross_beta*gamma_cross_beta;
	
/***----- Computation of the Energy -----***/
    if ( FLAG_orbit_spin == 'p' )
    {
      energy_square = ( numer1 - 2.0*sqrt( numer2 ) ) / denom1;
	    *energy = sqrt( energy_square );
    }
    else 
    {
      energy_square = ( numer1 + 2.0*sqrt( numer2 ) ) / denom1;	  
	    *energy = sqrt( energy_square );
    }

/**----- Computing the Angular Momentum component in the direction of the BH spin -----**/
    if ( FLAG_Schwarzschild == 'y' )
      angular_momentum_square = one_minus_z_m*( 3.0*rp2*energy_square - (*r_peri)*(3.0*(*r_peri)-4.0) );
	else
      angular_momentum_square = - one_minus_z_m*( alpha_cross_beta*energy_square + lambda_cross_beta ) / gamma_cross_beta;
  
    if ( FLAG_orbit_spin == 'p' )
      *Lz = sqrt( angular_momentum_square );
    else
      *Lz = - sqrt( angular_momentum_square );
	
/**----- Computing the Carter constant(s): C and Q -----**/
    if ( FLAG_Schwarzschild == 'y' )
      *q_carter = 0.0;
	else
      *q_carter = (*z_minus)*( -1.0*( alpha_cross_beta*energy_square + lambda_cross_beta)/gamma_cross_beta + spin2*( 1.0 - energy_square ) );
  
      *c_carter = *q_carter + ( spin*(*energy) - *Lz )*( spin*(*energy) - *Lz );


/**----- Computation of the Inclination Angle iota -----**/
    tani = sqrt(*c_carter)/(*Lz);
      
    cosi = *Lz/(sqrt(angular_momentum_square + *q_carter));
    
    *tilt = acos(cosi);
  }
/*----- Eccentric Case ( r_dot not zero ) -----*/
  else
  {   
/**----- Setting the values of the apocenter and pericenter (R_peri and R_apo) -----**/
    *r_peri = p_orbit/(1.0 + eccentricity);
    *r_apo  = p_orbit/(1.0 - eccentricity);
 
/**----- Constructing Auxiliary Quantities -----**/
    rp2 = (*r_peri)*(*r_peri);                               // r_peri^2
    ra2 = (*r_apo)*(*r_apo);                                 // r_apo^2
    one_minus_z_m = 1.0 - *z_minus;                          // 1 - z-
    two_mass_rp = 2.0*(*r_peri);                             // 2 M r_peri
    two_mass_ra = 2.0*(*r_apo);                              // 2 M r_apo
    delta_p = rp2 + spin2 - two_mass_rp;                     // r_peri^2 - 2 M r_peri + a^2
    delta_a = ra2 + spin2 - two_mass_ra;                     // r_apo^2 - 2 M r_apo + a^2
  
/**----- Constructing the coefficients of the equations for (E,Lz) [after substituting C(E,Lz)] -----**/
    alpha_a = (ra2 + spin2)*(ra2 + spin2*(*z_minus)) + two_mass_ra*spin2*one_minus_z_m;
    alpha_p = (rp2 + spin2)*(rp2 + spin2*(*z_minus)) + two_mass_rp*spin2*one_minus_z_m;
  
    beta_a = -two_mass_ra*spin;
    beta_p = -two_mass_rp*spin;
  
    gamma_a = - ( ra2 + spin2*(*z_minus) - two_mass_ra );
    gamma_p = - ( rp2 + spin2*(*z_minus) - two_mass_rp );
  
    lambda_a = - ( ra2 + spin2*(*z_minus) )*delta_a;
    lambda_p = - ( rp2 + spin2*(*z_minus) )*delta_p;
  
/**----- 'Cross Products': Coefficients of the combined equations: A_cross_B = A_a B_p - A_p B_a -----**/  
    alpha_cross_gamma = alpha_a*gamma_p - alpha_p*gamma_a;
    alpha_cross_beta = alpha_a*beta_p - alpha_p*beta_a;
    gamma_cross_beta = gamma_a*beta_p - gamma_p*beta_a;
    lambda_cross_gamma = lambda_a*gamma_p - lambda_p*gamma_a;
    lambda_cross_beta = lambda_a*beta_p - lambda_p*beta_a;
  
/**----- Computing the Energy -----**/
/***----- Auxiliary Quantities: E^2 = ( numer1 +- 2 sqrt(numer2) )/denom1 -----***/
    numer1 = - ( alpha_cross_gamma*lambda_cross_gamma + 2.0*one_minus_z_m*gamma_cross_beta*lambda_cross_beta );
  
    numer2 = one_minus_z_m*gamma_cross_beta*( lambda_cross_beta*( alpha_cross_gamma*lambda_cross_gamma 
	       + one_minus_z_m*gamma_cross_beta*lambda_cross_beta ) - alpha_cross_beta*lambda_cross_gamma*lambda_cross_gamma );
		 
    denom1 = alpha_cross_gamma*alpha_cross_gamma + 4.0*one_minus_z_m*alpha_cross_beta*gamma_cross_beta;

/***----- Computation of the Energy -----***/
    if ( FLAG_orbit_spin == 'p' )
    {
      energy_square = ( numer1 - 2.0*sqrt( numer2 ) ) / denom1;
	   *energy = sqrt( energy_square );
    }
    else 
    {
      energy_square = ( numer1 + 2.0*sqrt( numer2 ) ) / denom1;	  
	   *energy = sqrt( energy_square );
    }

/**----- Computing the Angular Momentum component in the direction of the BH spin -----**/
    if ( FLAG_Schwarzschild == 'y' )
      angular_momentum_square = one_minus_z_m*( (ra2 + (*r_peri)*(*r_apo) +rp2)*(energy_square-1.0) + 2.0*( *r_peri + *r_apo ) );
	else
      angular_momentum_square = - one_minus_z_m*( alpha_cross_beta*energy_square + lambda_cross_beta ) / gamma_cross_beta;
  
    if ( FLAG_orbit_spin == 'p' )
      *Lz = sqrt( angular_momentum_square );
    else 
      *Lz = - sqrt( angular_momentum_square );
	
/**----- Computing the Carter constant(s): C and Q -----**/
    if ( FLAG_Schwarzschild == 'y' )
	  *q_carter = 0.0; 
	else
      *q_carter = (*z_minus)*( -1.0*( alpha_cross_beta*energy_square + lambda_cross_beta )/gamma_cross_beta + spin2*( 1.0 - energy_square ) );
  
    *c_carter = *q_carter + ( spin*(*energy) - *Lz )*( spin*(*energy) - *Lz );

/**----- Computation of the Inclination Angle iota -----**/
    tani = sqrt(*q_carter)/(*Lz);
     
    cosi = *Lz/(sqrt(angular_momentum_square + *q_carter));
      
    *tilt = acos(cosi);
  }


/*----- This is the end -----*/
  return 0;
}

/*****************************************************************************************/
/*                                                                                       */
/*   The conventions of the ROK Code for the Carter Constant are:                        */
/*                                                                                       */
/*      4 .2       2   2           2                         2   2                       */
/*   rho  r  = [ (r + a )E - a Lz ]  - Delta [ Q + (a E - Lz) + r  ]                     */
/*                                                                                       */
/*                                                                                       */
/*   [there is an ambiguity between Q and C in this explanation]                         */
/*****************************************************************************************/
/*   NOTE: There are two 'Carter constants', namely Q and C.  They are                   */
/*   related by:                                                                         */
/*                                                                                       */
/*            Q = C + ( a E - Lz )^2                                                     */
/*                                                                                       */
/*   To see how they appear in the geodesic equations, have a look at the beginning of   */
/*   the file 'turning_points.c'.                                                        */
/*                                                                                       */
/*   In order to obtain (E,Lz,C/Q) from (p,e,theta_inc), this routine uses a procedure   */
/*   devised by the first time in a paper by W. Schmidt:                                 */
/*                                                                                       */
/*   W. Schmidt, "Celestial mechanics in Kerr spacetime", Classical and Quantum Gravity, */
/*   19, 2743–2764 (2002).                                                               */
/*                                                                                       */
/*                                                                                       */
/*---------------------------------------------------------------------------------------*/
/*                                                                                       */
/*   NOTE: The inverse function, that is, the one that goes from the constants of motion */
/*   (E,Lz,C/Q) to the orbital parameters (p,e,theta_inc) is contained in the c file:    */
/*                                                                                       */
/*                   ELQ_to_pez.c                                                        */
/*                                                                                       */
/*---------------------------------------------------------------------------------------*/
/*                                                                                       */
/*                         SOME TECHNICAL DETAILS                                        */
/*                                                                                       */
/*  We are going to deal first with the eccentric case (i.e., e != 0).  We will see      */
/*  later that the circular case is a singular case but that can be dealt with in the    */
/*  same way.                                                                            */
/*                                                                                       */
/* @ The equations for (E,Lz,C), which we call (eq_1,eq_2,eq_3), are obtained from the   */
/*   imposition that r_peri [ = p M / (1 + e)], r_apo [ = p M / (1 - e)], and z_-        */
/*   [ = cos^2(theta_min)] are extrema of the radial and polar motions respectively.   */
/*   Then, we can write:                                                                 */
/*                                                                                       */
/*     eq_1 = (d r/d tau)|r=r_peri = 0 ,                                                 */
/*                                                                                       */
/*     eq_2 = (d r/d tau)|r=r_apo = 0 ,                                                  */
/*                                                                                       */
/*     eq_3 = (d theta/ d tau)|theta=theta_min = 0 ,                                   */
/*                                                                                       */
/*   From the form of the Kerr geodesic equations in Boyer-Lindquist coordinates we can  */
/*   the following form for these equations:                                             */
/*                                                                                       */
/*     eq_1 = p1 E^2 + p2 E Lz + p3 Lz^2 + p4 C + p5 = 0 ,                               */
/*                                                                                       */
/*     eq_2 = q1 E^2 + q2 E Lz + q3 Lz^2 + q4 C + q5 = 0 ,                               */
/*                                                                                       */
/*     eq_3 = C - z_ [ Lz^2/(1-z_) + a^2 (1-E^2) ] = 0 ,                                 */
/*                                                                                       */
/*   where p1,p2,p3,p4,p5,q1,q2,q3,q4, and q5 are functions ONLY of the Black Hole       */
/*   parameters (M,a) and the turning points of the motion (r_peri, r_apo, z_-).  We can */
/*   see the system of equations (eq_1,eq_2,eq_3) as a system of equations for the       */
/*   unknowns (E,Lz,C).  One way to solve it is first to eliminate the unknown C, by     */
/*   isolating it from eq_3 and substituting it into equations eq_1 and eq_2.  Then, we  */
/*   we get two equations of the form:                                                   */
/*                                                                                       */
/*    eq_4 = alpha_a E^2 + 2 beta_a E Lz + gamma_a Lz^2 + lambda_a = 0 ,                 */
/*                                                                                       */
/*    eq_5 = alpha_p E^2 + 2 beta_p E Lz + gamma_p Lz^2 + lambda_p = 0 ,                 */
/*                                                                                       */
/*   where the coefficients alpha_I, beta_I, gamma_I, and lambda_I (I = a,p; where 'a'   */
/*   stands for apocenter and 'p' for pericenter) are given by                           */
/*                                                                                       */
/*   alpha_I = ( r_I^2 + a^2 ) ( r_I^2 + a^2 z_- ) + 2 M r_I a^2 ( 1 - z_- ) ,           */
/*                                                                                       */
/*   beta_I = - 2 M r_I a ,                                                              */
/*                                                                                       */
/*   gamma_I = - [ r_I^2 + a^2 z_- - 2 M r_I ] / ( 1 - z_- ) ,                           */
/*                                                                                       */
/*   lambda_I - ( r_I^2 + a^2 ) ( r_I^2 - 2 M r_I + a^2 ) ,                              */
/*                                                                                       */
/*   and as it is clear from the notation, when I = a, means to evaluate these           */
/*   coefficients for r = r_apo, and when I = p, means to evaluate them for r = r_peri.  */
/*   At this point, we can combine equations (eq_4,eq_5) in order to eliminate one of    */
/*   the two unknowns, for instance Lz.  The resulting equation for E, the energy, is a  */
/*   biquadratic equation of the form:                                                   */
/*                                                                                       */
/*   ( |alpha,gamma|^2 + 4 |alpha,beta| |gamma,beta| ) E^4                               */
/*   + 2 ( |alpha,gamma| |lambda,gamma| + 2 |gamma,beta| |lambda,beta| ) E^2             */
/*   + |lambda,gamma|^2 = 0 ,                                                            */
/*                                                                                       */
/*   where we have used the following notation                                           */
/*                                                                                       */
/*      |A,B| := A_a B_p - A_p B_a ,                                                     */
/*                                                                                       */
/*   where again the subscripts 'a' and 'p' mean evaluation at r = r_apo and r = r_peri  */
/*   respectively.  We can see this equation a normal quadratic equation for E^2. Then,  */
/*   it has two solutions for E^2, that can be written as:                               */
/*                                                                                       */
/*                            - C_1 +- 2 Sqrt[ C_2 ]                                     */
/*                     E^2 = ------------------------- ,                                 */
/*                                     C_3                                               */
/*                                                                                       */
/*   where the quantities C_1, C_2, and C_3 are given by                                 */
/*                                                                                       */
/*   C_1 = |alpha,gamma| |lambda,gamma| + 2 |gamma,beta| |lambda,beta| ,                 */
/*                                                                                       */
/*   C_2 = |gamma,beta| { |lambda,beta| ( |alpha,gamma| |lambda,gamma|                   */
/*       + |gamma,beta| |lambda,beta| ) - |alpha,beta| |lambda,gamma|^2 } ,              */
/*                                                                                       */
/*   C_3 = |alpha,gamma|^2 + 4 |alpha,beta| |gamma,beta| .                               */
/*                                                                                       */
/*   The sign '+' corresponds to retrograde orbits as these ones have the greatest       */
/*   energy, whereas the sign '-' corresponds to prograde orbits.  Given one of these    */
/*   two solutions or E^2 we still have the sign of Sqrt[ E^2 ]. We take the positive    */
/*   sign, as the negative one just represents the time reverse motion.   Then, once we  */
/*   have solve for E, we can find the angular momentum through the following expression */
/*                                                                                       */
/*       Lz^2 = { |alpha,beta| E^2 + |lambda,beta| } / |beta,gamma| .                    */
/*                                                                                       */
/*   Here we have also two possible signs for Lz.  It is simple to see that the positive */
/*   sign corresponds to prograde orbits whereas the negative one corresponds to         */
/*   retrograde orbits.    With E and Lz, we can compute C from eq_3 above, that is,     */
/*                                                                                       */
/*                                                                                       */
/*                   Lz^2        2        2                                              */
/*       Q = z_- { ---------- + a  ( 1 - E  ) }   .                                      */
/*                  1 - z_-                                                              */
/*                                                                                       */
/*   and obviously we can compute Q by usign the definition                              */
/*                                                                                       */
/*       Q = C - ( a E - Lz )^2 ,                                                        */
/*                                                                                       */
/*   and the tilt angle iota by                                                   */
/*                                                                                       */
/*                           Lz                                                          */
/*       cos(iota) = -------------------- .                                              */
/*                     Sqrt[ Lz^2 + C ]                                                  */
/*                                                                                       */
/*                                                                                       */
/*   Finally, one can find the other extremal points (have a look at the C function      */
/*   'compute_extrema_of_motion.c') of the radial and polar motion (r_3, r_4, z_+) from  */
/*   the equations that they satisfy.  For radial motion these equations are             */
/*                                                                                       */
/*     r_peri + r_apo + r_3 + r_4 = 2 M / ( 1 - E^2 ) ,                                  */
/*                                                                                       */
/*     r_peri r_apo + r_3 r_4 + ( r_3 + r_4 ) ( r_peri + r_apo )                         */
/*       = [ a^2 ( 1 - E^2 ) + Lz^2 + C ] / ( 1 - E^2 ) ,                                */
/*                                                                                       */
/*     r_3 r_4 ( r_peri + r_apo ) + r_peri r_apo ( r_3 + r_4 ) = 2 M Q / ( 1 - E^2 ) ,   */
/*                                                                                       */
/*     r_peri r_apo r_3 r_4 = a^2 C / ( 1 - E^2 ) ,                                      */
/*                                                                                       */
/*   and the equations for the case of the polar motion are                              */
/*                                                                                       */
/*     z_- + z_+ = [ C + Lz^2 + a^2 ( 1 - E^2 ) ] / [ a^2 ( 1 - E^2 ) ] ,                */
/*                                                                                       */
/*     z_- z_+ = C / [ a^2 ( 1 - E^2 ) ] .                                               */
/*                                                                                       */
/*  Since (r_peri,r_apo) are known, we need two use just two of the four equations for   */
/*  the extrema of the radial motion to get (r_3,r_4).  We have used the first one and   */
/*  the fourth one.  Then, one case use (as we do in this code) the other two to check   */
/*  the results.  Regarding the polar motion, we just need to use one of the equations   */
/*  for the extrema.  We use the first one, then we can use (as we do in this code) the  */
/*  second one to check the results.                                                     */
/*                                                                                       */
/*  Finally, let us comment on the case of circular (spherical) orbits.  This case is    */
/*  characterized by the vanishing of the eccentricity, which means that r_peri = r_apo. */
/*  Then, we have: alpha_a = alpha_p, beta_a = beta_p, gamma_a = gamma_p, and lambda_a = */
/*  lambda_p.  That is, the two equations for (E,Lz) are the same and we cannot isolate  */
/*  neither E nor Lz.  The way out comes from the characterization of the circular       */
/*  orbits in terms of the function R(r):                                                */
/*                                                                                       */
/*    R(r_o) = 0,  R'(r_o) = 0, and R''(r_o) > 0.   [ r_o = r_peri = r_apo ]             */
/*                                                                                       */
/*  The first equation coincide with the two equations of the eccentric case (which are  */
/*  the same now as we have already mentioned).  The second one is the extra equation    */
/*  that we need in order to have two equations for two unknowns, (E,Lz).  Then, the     */
/*  coefficients alpha, beta, gamma, and lambda for this second equation of the circular */
/*  case are:                                                                            */
/*                                                                                       */
/*    alpha = 2 r_o ( r_o^2 + a^2 ) - a^2 ( r_o - M ) ( 1 - z_- ) ,                      */
/*                                                                                       */
/*    beta = - a M ,                                                                     */
/*                                                                                       */
/*    gamma = - ( r - M ) / ( 1 - z_- ) ,                                                */
/*                                                                                       */
/*    lambda = - r_o ( r_o^2 - 2 M r_o + a^2 ) - ( r_o - M ) ( r_o^2 + a^2 z_- ) .       */
/*                                                                                       */
/*  Then, following exactly the same procedure as described for the eccentric case we    */
/*  will obtain the constants of motion and the other relevant quantities.               */
/*                                                                                       */
/*****************************************************************************************/

