
/*==============================================================================
 Computing the Christoffel symbols of the Kerr metric in Boyer-Lindquist
 coordinates
 -------------------------------------------------------------------------------
 Created by Priscilla on 29/12/2012
 -------------------------------------------------------------------------------
 Last Update : 23.09.13
 =============================================================================*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "global_quantities_kerr.h"


int christoffel_symbols(double r, double costh, double gamma[4][4][4])
{
  double r2, r3, r4, r5,r2pa2;
  double f1, f2, f3, fr;
  double cos_theta, sin_theta;
  double costh2, sinth2, sincosth;
  
  double rho2, rho4, rho6;
  double Spin2, Spin3, Spin4;
  
/*----- Computing some Frequent Quantities -----*/
  r2 = r*r;
  r3 = r2*r;
  r4 = r2*r2;
  r5 = r4*r;

  Spin2 = SPIN*SPIN;
  Spin3 = Spin2*SPIN;
  Spin4 = Spin2*Spin2;
  
  r2pa2 = r2 + Spin2;
  
  cos_theta = costh;
  sin_theta = sqrt(1. - costh*costh);
  
  costh2 = cos_theta*cos_theta;
  sinth2 = sin_theta*sin_theta;
  sincosth = sin_theta*cos_theta;
  
  rho2 = r2 + Spin2*costh2;
  rho4 = rho2*rho2;
  rho6 = rho4*rho2;
    
  f1 = r2pa2 - 2.0*r;
  f2 = r2 - Spin2*costh2;
  f3 = 1.0/r;
  fr = 1.0 - 2.0*f3;
  
/*----- Computing the Different Christoffel Symbols -----*/
/*----- Notation: z = theta, f = phi                -----*/

/**----- G^t_tt -----**/
  gamma[0][0][0] = 0.0;               // gBL_ttt = 0 
  
/**----- G^t_tr -----**/
  gamma[0][0][1] = r2pa2*f2/(f1*rho4);
  
/**----- G^t_tz -----**/   
  gamma[0][0][2] = -2.0*Spin2*r*sincosth/rho4;
  
/**----- G^t_tf -----**/  
  gamma[0][0][3] = 0.0;               // gBL_ttf = 0
  
/**----- G^t_rr -----**/  
  gamma[0][1][1] = 0.0;               // gBL_trr = 0
  
/**----- G^t_rz -----**/  
  gamma[0][1][2] = 0.0;               // gBL_trz = 0
  
/**----- G^t_rf -----**/  
  gamma[0][1][3] = sinth2*( SPIN*( Spin2*(Spin2 - r2)*costh2  - r2*(Spin2+3.0*r2) )/(f1*rho4));
                 
/**----- G^t_zz -----**/  
  gamma[0][2][2] = 0.0;               // gBL_tzz = 0

/**----- G^t_zf -----**/  
  gamma[0][2][3] = 2.0*sincosth*sinth2*r*Spin3/rho4;  

/**----- G^t_ff -----**/  
  gamma[0][3][3] = 0.0;               // gBL_tff = 0

/**----- G^r_tt -----**/
  gamma[1][0][0] = f1*f2/rho6;
  
/**----- G^r_tr -----**/  
  gamma[1][0][1] = 0.0;               // gBL_rtr = 0

/**----- G^r_tz -----**/
  gamma[1][0][2] = 0.0;               // gBL_rtz = 0
  
/**----- G^r_tf -----**/  
  gamma[1][0][3] = sinth2*( -gamma[1][0][0]*SPIN  );

/**----- G^r_rr -----**/  
  gamma[1][1][1] = -( r*(r-Spin2) + Spin2*costh2*(r - 1.0) )/(f1*rho2);

/**----- G^r_rz -----**/
  gamma[1][1][2] = -Spin2*sincosth/rho2;

/**----- G^r_rf -----**/  
  gamma[1][1][3] = 0.0;               // gBL_rrf = 0 
  
/**----- G^r_zz -----**/  
  gamma[1][2][2] = -f1*r/rho2;

/**----- G^r_zf -----**/  
  gamma[1][2][3] = 0.0;               // gBL_rzf = 0
  
/**----- G^r_ff -----**/  
  gamma[1][3][3] = -f1*( Spin2*sinth2*(Spin2*costh2-r2) + r*rho4 )*sinth2/rho6;  

/**----- G^z_tt -----**/
  gamma[2][0][0] = -2.0*Spin2*r*sincosth/rho6;

/**----- G^z_tr -----**/  
  gamma[2][0][1] = 0.0;               // gBL_ztr = 0

/**----- G^z_tz -----**/  
  gamma[2][0][2] = 0.0;               // gBL_ztz = 0

/**----- G^z_tf -----**/  
  gamma[2][0][3] = 2.0*r*SPIN*r2pa2*sincosth/rho6;

/**----- G^z_rr -----**/  
  gamma[2][1][1] = Spin2*sincosth/(f1*rho2);

/**----- G^z_rz -----**/  
  gamma[2][1][2] = r/rho2;
  
/**----- G^z_rf -----**/  
  gamma[2][1][3] = 0.0;               // gBL_zrf = 0
  
/**----- G^z_zz -----**/  
  gamma[2][2][2] = -Spin2*sincosth/rho2;
  
/**----- G^z_zf -----**/  
  gamma[2][2][3] = 0.0;               // gBL_zzf = 0
  
/**----- G^z_ff -----**/  
  gamma[2][3][3] = -( 2.0*Spin2*r*sinth2*(2.0*r2+Spin2*(1+costh2)) + (Spin2+r2)*rho4 )*cos_theta*sin_theta/rho6;

/**----- G^f_tt -----**/
  gamma[3][0][0] = 0.0;               // gBL_ftt = 0
  
/**----- G^f_tr -----**/  
 gamma[3][0][1] = SPIN*f2/(f1*rho4);

/**----- G^f_tz -----**/
  gamma[3][0][2] = (cos_theta/sin_theta)*( -2.0*r*SPIN/rho4  );

/**----- G^f_tf -----**/  
  gamma[3][0][3] = 0.0;               // gBL_ftf = 0
  
/**----- G^f_rr -----**/  
  gamma[3][1][1] = 0.0;               // gBL_frr = 0
  
/**----- G^f_rz -----**/  
  gamma[3][1][2] = 0.0;               // gBL_frz = 0

/**----- G^f_rf -----**/  
  gamma[3][1][3] = ( r*rho4 + Spin4*costh2*sinth2 - r2*Spin2*(1.0+costh2) - 2.0*r4 )/(f1*rho4);
  
/**----- G^f_zz -----**/  
  gamma[3][2][2] = 0.0;               // gBL_fzz = 0

/**----- G^f_zf -----**/  
  gamma[3][2][3] = cos_theta*( 2.0*r*Spin2*sinth2 + rho4 )/(rho4*sin_theta);
  
/**----- G^f_ff -----**/  
  gamma[3][3][3] = 0.0;               // gBL_fff = 0


/*----- This is the end -----*/
  return 0;

}




