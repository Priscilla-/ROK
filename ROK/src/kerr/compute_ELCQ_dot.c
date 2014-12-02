/*==============================================================================
     Given the orbital parameters of a geodesic orbit around a Kerr BH,                  

                         (p,e,iota/theta_inc),

   which in particular means to know already (E,Lz,C/Q), this routine computes

                         (E_dot,L_dot,Q_dot),
 -------------------------------------------------------------------------------
                                            Created by Priscilla on 29/12/2012
 -------------------------------------------------------------------------------
  Last Update : 28.12.12
 =============================================================================*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "global_quantities_kerr.h"
#include "coefficients_of_the_fit.h"
#include "global_quantities_main.h"

#include "macros_2.h"

int fluxes_scalar_field(void);


int compute_ELCQ_dot(double p_orbit, double eccentricity, double inclination,double angular_momentum_z, double q_carter,
                     double energy_0e, double angular_momentum_z_0e, double q_carter_0e,
					 double *energy_dot, double *angular_momentum_z_dot, double *q_carter_dot)
{
  double eccentricity2,eccentricity4,eccentricity6;           // e^2, e^4, e^6
  double Spin3, Spin4;                                        // (a/M)^3, (a/M)^4
  double one_minus_ecc2;                                      // 1 - e^2
  double one_minus_ecc2_3halfs;                               // ( 1 - e^2 )^(3/2)
  double one_over_one_minus_ecc2_3halfs;                      // ( 1 - e^2 )^(-3/2)

  double p_2,p_4;                                             // p^2, p^4
  double one_over_p;                                          // 1 / p
  double one_over_p_sqrt;                                     // sqrt( 1 / p )
  double one_over_p_3halfs;                                   // ( 1 / p )^(3/2)
  double one_over_p_2;                                        // ( 1 / p )^2
  double one_over_p_5halfs;                                   // ( 1 / p )^(5/2)
  double one_over_p_7halfs;                                   // ( 1 / p )^(7/2)
  double one_over_p_4;                                        // ( 1 / p )^4
  double one_over_p_5;                                        // ( 1 / p )^5
  
  double cos_i,cos_i_2,cos_i_3,cos_i_4,cos_i_5;               // cos(Inclination), cos(Inclination)^2, cos(Inclination)^3, cos(Inclination)^4, cos(Inclination)^5
  double sin_i,sin_i_2;                                       // sin(Inclination), sin(Inclination)^2
  double tan_i;                                               // tan(Inclination)
  
  double g_1,g_2,g_3,g_4,g_5,g_6,g_9,ga_10,gb_10;             // g_i(e)
  double g_11,g_12,g_13,g_14;                                 // g_i(e)
  double g0_1,g0_2,g0_3,g0_4,g0_5,g0_6,g0_9,ga0_10,gb0_10;    // g_i(0)
  double g0_11,g0_12,g0_13,g0_14;                             // g_i(0)
  
  double n1_0;                                                // N1(p,i) = N1(r_apo(e=0)) = E(e=0,p,i) ( r0^4 + a^2 r0^2 ) - 2 a M [ Lz(e=0,p,i) - a E(e=0,p,i) ] r0
  double n4_0;                                                // N4(p,i) = M^2 ( 2 p - p^2 ) Lz(e=0,p,i) - 2 M^2 a p E(e=0,p,i)
  double n5_0;                                                // N5(p,i) = M^2 ( 2 p - p^2 - (a/M)^2 ) / 2  
  
  double e_dot_2PN;                                           // 2 PN Approximation to the Energy Flux
  double e_dot_2PN_0e;                                        // 2 PN Approximation to the Energy Flux in the case e = 0
  double lz_dot_2PN;                                          // 2 PN Approximation to the Angular Momentum (in the Spin Direction) Flux
  double lz_dot_2PN_0e;                                       // 2 PN Approximation to the Angular Momentum (in the Spin Direction) Flux
  double lz_dot_fit;                                          // Fit of the Angular Momentum (in the Spin Direction) Flux to the corresponding Teukolsky Flux
  double q_dot_over_sqrtQ_2PN;                                // 2 PN Approximation to the Carter Constant Flux
  double q_dot_over_sqrtQ_2PN_0e;                             // 2 PN Approximation to the Carter Constant Flux in the case e = 0
  double inclination_dot_fit;                                 // Fit of the 'Inclination Flux' to the corresponding Teukolsky Flux
  
  double factor_sin_i_2_over_sqrt_Q;
  double angular_momentum_z_factor;
  double inclination_dot_fit_without_factor;
    

/*----- 2PN Expressions -----*/
/**----- Some coefficients/quantities -----**/
  eccentricity2 = SQR(eccentricity);
  eccentricity4 = SQR(eccentricity2);
  eccentricity6 = SQR(eccentricity2);

  Spin3 = SPIN2*SPIN;
  Spin4 = SQR(SPIN2);
    
  one_minus_ecc2 = 1.0 - eccentricity2;                       
  one_minus_ecc2_3halfs = pow(one_minus_ecc2,1.5);
  one_over_one_minus_ecc2_3halfs = 1.0/one_minus_ecc2_3halfs;
  
  p_2 = SQR(p_orbit);
  p_4 = SQR(p_2);
  one_over_p = 1.0/p_orbit;
  one_over_p_sqrt = sqrt(one_over_p);
  one_over_p_3halfs = one_over_p_sqrt*one_over_p;
  one_over_p_2 = SQR(one_over_p);
  one_over_p_5halfs = one_over_p_3halfs*one_over_p;
  one_over_p_7halfs = one_over_p_3halfs*one_over_p_2;
  one_over_p_4 = SQR(one_over_p_2);
  one_over_p_5 = one_over_p_4*one_over_p;
  
  cos_i = cos(inclination);
  cos_i_2 = SQR(cos_i);
  cos_i_3 = cos_i_2*cos_i;
  cos_i_4 = SQR(cos_i_2);
  cos_i_5 = cos_i_3*cos_i_2;
  sin_i = sin(inclination);
  sin_i_2 = SQR(sin_i);
  tan_i = sin_i/cos_i;                                     
 
/**----- Functions evaluated for circular orbits -----**/
  n1_0 = energy_0e*( p_4 + SPIN2*p_2 ) - 2.0*SPIN*p_orbit*( angular_momentum_z_0e - SPIN*energy_0e );  
  n4_0 = ( 2.0*p_orbit - p_2 )*angular_momentum_z_0e - 2.0*SPIN*p_orbit*energy_0e; 
  n5_0 = 0.5*( 2.0*p_orbit - p_2 - SPIN2 ); 
      
/**----- g_i(e) functions -----**/
  g_1  = 1.0 + (73.0/24.0)*eccentricity2  + (37.0/96.0)*eccentricity4;
  g_2  = 73.0/12.0 + (823.0/24.0)*eccentricity2 + (949.0/32.0)*eccentricity4 + (491.0/192.0)*eccentricity6;  
  g_3  = 1247.0/336.0 + (9181.0/672.0)*eccentricity2;
  g_4  = 4.0 + (1375.0/48.0)*eccentricity2;
  g_5  = 44711.0/9072.0 + (172157.0/2592.0)*eccentricity2;
  g_6  = 33.0/16.0 + (359.0/32.0)*eccentricity2;
  g_9  = 1.0 + (7.0/8.0)*eccentricity2;
  ga_10 = 61.0/24.0 + (63.0/8.0)*eccentricity2 + (95.0/64.0)*eccentricity4;
  gb_10 = 61.0/8.0 + (91.0/4.0)*eccentricity2 + (461.0/64.0)*eccentricity4;
  g_11 = 1247.0/336.0 + (425.0/336.0)*eccentricity2;
  g_12 = 4.0 + (97.0/8.0)*eccentricity2;
  g_13 = 44711.0/9072.0 + (302893.0/6048.0)*eccentricity2;
  g_14 = 33.0/16.0 + (95.0/16.0)*eccentricity2;

/**----- g_i(0) functions -----**/
  g0_1  = 1.0;
  g0_2  = 73.0/12.0;
  g0_3  = 1247.0/336.0;
  g0_4  = 4.0;
  g0_5  = 44711.0/9072.0;
  g0_6  = 33.0/16.0;
  g0_9  = 1.0;
  ga0_10 = 61.0/24.0;
  gb0_10 = 61.0/8.0;
  g0_11 = g0_3;
  g0_12 = 4.0;
  g0_13 = g0_5;
  g0_14 = g0_6;

/**----- Fluxes -----**/
  e_dot_2PN = -(32.0/5.0)*MASS_ratio2*one_over_p_5*one_minus_ecc2_3halfs*( 
			 g_1 - SPIN*one_over_p_3halfs*g_2*cos_i - one_over_p*g_3 
		   + PI*one_over_p_3halfs*g_4 - one_over_p_2*g_5 + SPIN2*one_over_p_2*g_6
           - ( (527.0/96.0) + (6533.0/192.0)*eccentricity2 )*SPIN2*one_over_p_2*sin_i_2 );

  lz_dot_2PN = -(32.0/5.0)*MASS_ratio2*one_over_p_7halfs*one_minus_ecc2_3halfs*(  
              g_9*cos_i + SPIN*one_over_p_3halfs*(ga_10 - cos_i_2*gb_10) 
			- one_over_p*g_11*cos_i + PI*one_over_p_3halfs*g_12*cos_i
			- one_over_p_2*g_13*cos_i + SPIN2*one_over_p_2*cos_i*(g_14 - ( (45.0/8.0) + 18.5*eccentricity2 )*sin_i_2) );
  
  q_dot_over_sqrtQ_2PN = -(64.0/5.0)*MASS_ratio2*one_over_p_7halfs*sin_i*one_minus_ecc2_3halfs*(
                         g_9 - SPIN*one_over_p_3halfs*cos_i*gb_10 - one_over_p*g_11
                       + PI*one_over_p_3halfs*g_12 - one_over_p_2*g_13
                       + SPIN2*one_over_p_2*(g_14 - ( (45.0/8.0) )*sin_i_2) );

/**----- Fluxes at e=0 -----**/
  e_dot_2PN_0e = -(32.0/5.0)*MASS_ratio2*one_over_p_5*( g0_1 - SPIN*one_over_p_3halfs*g0_2*cos_i 
               - one_over_p*g0_3 + PI*one_over_p_3halfs*g0_4 - one_over_p_2*g0_5 + SPIN2*one_over_p_2*g0_6
               - (527.0/96.0)*SPIN2*one_over_p_2*sin_i_2 );

  lz_dot_2PN_0e = -(32.0/5.0)*MASS_ratio2*one_over_p_7halfs*(  
                  g0_9*cos_i + SPIN*one_over_p_3halfs*(ga0_10 - cos_i_2*gb0_10) 
			    - one_over_p*g0_11*cos_i + PI*one_over_p_3halfs*g0_12*cos_i
			    - one_over_p_2*g0_13*cos_i + SPIN2*one_over_p_2*cos_i*(g0_14 - (45.0/8.0)*sin_i_2) );
  
  q_dot_over_sqrtQ_2PN_0e = -(64.0/5.0)*MASS_ratio2*one_over_p_7halfs*sin_i*(
                            g0_9 - SPIN*one_over_p_3halfs*cos_i*gb0_10 - one_over_p*g0_11
                          + PI*one_over_p_3halfs*g0_12 - one_over_p_2*g0_13
                          + SPIN2*one_over_p_2*(g0_14 - (45.0/8.0)*sin_i_2) );


/*----- Expressions fitted to Teukolsky Fluxes -----*/
/**----- Fluxes -----**/
  lz_dot_fit = -(32.0/5.0)*MASS_ratio2*one_over_p_7halfs*( 
               cos_i + SPIN*one_over_p_3halfs*( 61.0/24.0 - (61.0/8.0)*cos_i_2 ) - (1247.0/336.0)*one_over_p*cos_i
			 + 4.0*PI*one_over_p_3halfs*cos_i - (44711.0/9072.0)*one_over_p_2*cos_i
			 + SPIN2*one_over_p_2*cos_i*(33.0/16.0 - (45.0/8.0)*sin_i_2) 
			 + one_over_p_5halfs*( SPIN*(DA_1 + DB_1*one_over_p_sqrt + DC_1*one_over_p) + Spin3*(DA_2 + DB_2*one_over_p_sqrt + DC_2*one_over_p)
			                     + cos_i*(CA_1 + CB_1*one_over_p_sqrt + CC_1*one_over_p) + SPIN2*cos_i*(CA_2 + CB_2*one_over_p_sqrt + CC_2*one_over_p)
                                 + Spin4*cos_i*(CA_3 + CB_3*one_over_p_sqrt + CC_3*one_over_p) + SPIN*cos_i_2*(CA_4 + CB_4*one_over_p_sqrt + CC_4*one_over_p)
								 + Spin3*cos_i_2*(CA_5 + CB_5*one_over_p_sqrt + CC_5*one_over_p) + SPIN2*cos_i_3*(CA_6 + CB_6*one_over_p_sqrt + CC_6*one_over_p)
                                 + Spin4*cos_i_3*(CA_7 + CB_7*one_over_p_sqrt + CC_7*one_over_p) + Spin3*cos_i_4*(CA_8 + CB_8*one_over_p_sqrt + CC_8*one_over_p)
                                 + Spin4*cos_i_5*(CA_9 + CB_9*one_over_p_sqrt + CC_9*one_over_p) ) 
			 + one_over_p_7halfs*SPIN*cos_i*( FA_1 + FB_1*one_over_p_sqrt + SPIN*(FA_2 + FB_2*one_over_p_sqrt) + SPIN2*(FA_3 + FB_3*one_over_p_sqrt)
			                                + cos_i_2*(FA_4 + FB_4*one_over_p_sqrt) + SPIN*cos_i_2*(FA_5 + FB_5*one_over_p_sqrt) 
											+ SPIN2*cos_i_2*(FA_6 + FB_6*one_over_p_sqrt) ) );
  if ( sin_i_2 < 1.0e-15 )
  {
    factor_sin_i_2_over_sqrt_Q = fabs(sin_i*cos_i/angular_momentum_z);
  }
  else
  {
    factor_sin_i_2_over_sqrt_Q = sin_i_2/sqrt(q_carter);
  }

  inclination_dot_fit_without_factor = (32.0/5.0)*MASS_ratio2*SPIN*one_over_p_5*( 61.0/24.0 
                                     + one_over_p*(DA_1 + DB_1*one_over_p_sqrt + DC_1*one_over_p) + SPIN2*one_over_p*(DA_2 + DB_2*one_over_p_sqrt + DC_2*one_over_p)  
					                 + SPIN*cos_i*one_over_p_sqrt*(CA_10 + CB_10*one_over_p + CC_10*one_over_p_3halfs) 
                                     + SPIN2*cos_i_2*one_over_p*(CA_11 + CB_11*one_over_p_sqrt + CC_11*one_over_p)  
					                 + one_over_p_5halfs*Spin3*cos_i*( FA_7 + FB_7*one_over_p_sqrt + SPIN*(FA_8 + FB_8*one_over_p_sqrt) 
													                 + cos_i_2*(FA_9 + FB_9*one_over_p_sqrt) + SPIN*cos_i_2*(FA_10 + FB_10*one_over_p_sqrt) ) );

  inclination_dot_fit = factor_sin_i_2_over_sqrt_Q*inclination_dot_fit_without_factor;
  		 											  

/*----- FINAL EXPRESSIONS FOR THE FLUXES -----*/	
    
  *energy_dot = one_minus_ecc2_3halfs*( one_over_one_minus_ecc2_3halfs*e_dot_2PN
                - e_dot_2PN_0e - (n4_0/n1_0)*lz_dot_2PN_0e
                - (n5_0/n1_0)*q_dot_over_sqrtQ_2PN_0e*sqrt(q_carter_0e) )/MASS_ratio;
							  
  *angular_momentum_z_dot = one_minus_ecc2_3halfs*( one_over_one_minus_ecc2_3halfs*lz_dot_2PN
                           - lz_dot_2PN_0e + lz_dot_fit )/MASS_ratio;

  if ( sin_i_2 < 1.0e-15 )
  {
     angular_momentum_z_factor = angular_momentum_z_0e/angular_momentum_z;

     *q_carter_dot = one_minus_ecc2_3halfs*sqrt(q_carter)*( one_over_one_minus_ecc2_3halfs*q_dot_over_sqrtQ_2PN
                    - q_dot_over_sqrtQ_2PN_0e
                     + 2.0*tan_i*( lz_dot_fit + angular_momentum_z_factor*inclination_dot_fit_without_factor ) )/MASS_ratio;
  }
  else
  {
     *q_carter_dot = one_minus_ecc2_3halfs*sqrt(q_carter)*( one_over_one_minus_ecc2_3halfs*q_dot_over_sqrtQ_2PN
                    - q_dot_over_sqrtQ_2PN_0e
                    + 2.0*tan_i*( lz_dot_fit + sqrt(q_carter_0e)*inclination_dot_fit/sin_i_2 ) )/MASS_ratio;
  }
 
/*----- This is the end -----*/
  return 0;

}
/*****************************************************************************************/
/*                                                                                       */
/*  NOTE:  The conventions of the Code for the Carter Constant are:                      */
/*                                                                                       */
/*      4 .2       2   2           2                         2   2                       */
/*   rho  r  = [ (r + a )E - a Lz ]  - Delta [ Q + (a E - Lz) + r  ]                     */
/*                                                                                       */
/*   Q = C + ( a E - Lz )^2                                                              */
/*                                                                                       */
/*****************************************************************************************/


