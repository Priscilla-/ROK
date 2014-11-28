/*=====================================================================================
 
   This routine computes analytical solutions of bound timelike geodesics in  Kerr in
   function of the "Mino" time lambda for MBH spin smaller than MBH mass.
 
   The expressions computed here are quoted in Class. Quantum Grav. 26 (2009) 135002 [1], 
   where the radial solution r(lamba) is obtained by inverting lambda(r) from the geode-
   sic equations and theta(lambda) is obtained by inverting lambda(theta). 
   
   To perform their computation we need to specify r(0) and theta(0) as well as the sign
   of their first derivatives with respect lambda. Here we set r(0)= r_peri = pM/(1+e) and 
   theta(0) = pi/2, so we can use the solutions quoted in [1] independently of the sign 
   of their derivatives, and the initial angle variables are psi = 0 = xi.

 --------------------------------------------------------------------------------------
                                                    Created by Priscilla on 20/12/2012
 --------------------------------------------------------------------------------------
   Last Update : 27.12.12
 ======================================================================================*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "global_quantities_main.h"
#include "physical_quantities.h"
#include "global_quantities_kerr.h"

#include "macros_1.h"
#include "macros_2.h"

double ellK(double, double);

void JacobianEllipticFunctions(double, double, double *, double * , double * );


int analytical_trajectory(int l_time, double dMt)
{
  FILE *data;
  double sn_ur, sn_uth;                     // Jacobi's elliptic function and its square for the radial and polar solutions
  double r_M, theta_M;                      // Radial and Polar coordinates in "Mino" time
  double ur_M, uth_M;                       // Arguments for the Jacobi's elliptic function
  double Kcomp_r, Kcomp_th;                 // Complete elliptic integral of the firs kind
  double kr, kr_num, kr_dem;                // Arguments of the complete elliptic integrals
  double kth;                               // Argument of the complete elliptic integrals
  double yr, yth;                           // Arguments of the incomplete ellictic integral of first kind: arcsin(phy)
  double sn, cn, dn;                        // Arguments of Jacobi's elliptic function
  double Fin_r, Fin_th;                     // Incomplete elliptical functions of the first kind for the radial an polar solutions
  double lambda1_r, lambda1_th;             // Auxiliar times
  double epsilon0, M_r1, M_r2;              // Auxiliar variables
  double M_th1, M_th2, dMt_r, dMt_th;       // Auxiliar variables
  double omega_r, omega_th;                 // Auxiliar variables
  double lambda_r, lambda_th;                 // Auxiliar variables
  double mode_th, mode_r;                   // Auxiliar variables      
    
   
    
    double rho2;                                   // rho^2 = r^2 + a^2cos^2(theta)
    double delta;                                  // delta = r^2 - 2*Mr + a^2
    double dT_p, dTau;                             // BL time and proper time: dTau = rho^2*dT_p/T(r,Theta)
    double sigma2;                                 // sigma^2 = (r^2 +a^2)^2 -a^2*Delta*sin^2(theta)
    double dMt;                                    // Resonance time
    // int l_time;                                    // Counter
    
/*----Computing fundamental frequencies ---*/
    if (GV_Fundamental_Frequencies =='y')
    {
        rho2 = SQR(GV_R_p) + GV_Spin2*SQR(cos(GV_Theta_p));
        delta = SQR(GV_R_p) - 2.0*GV_R_p + GV_Spin2;
        sigma2 = SQR( SQR(GV_R_p) + GV_Spin2)- GV_Spin2*delta*SQR(sin(GV_Theta_p));
        
        dT_p = n_global*GV_Dt_ini;
        dTau = dT_p*rho2*delta/(sigma2*GV_E - 2.0*GV_Spin*GV_R_p*GV_Lz);
        
        dMt = dTau/rho2;
        
        // compute_frequencies();
        
        //analytical_trajectory(l_time, dMt);
        
        l_time ++;
    }

    
    
  omega_r = GV_Omega_r;
  omega_th = GV_Omega_th;
  
  lambda_r = 2.0*PI/omega_r;
  lambda_th = 2.0*PI/omega_th;
    
/*---Computing the radial and polar times---*/
  dMt_r =  dMt - floor(dMt/lambda_r)*lambda_r;
  dMt_th = dMt - floor(dMt/lambda_th)*lambda_th;
    
/*---Computing arguments of the complete elliptical integrals---*/
  kr_num = (GV_R_apo - GV_R_peri)*(GV_R_3 - GV_R_4);
  kr_dem = (GV_R_apo - GV_R_3)*(GV_R_peri - GV_R_4);
    
  mode_r = kr_num/kr_dem;
  kr = sqrt(mode_r);
    
  mode_th = GV_Z_minus/GV_Z_plus;
  kth = sqrt(mode_th);
    
/*---Computing Legendre elliptical integrals First kind---*/
  Kcomp_r  = ellK(0.5*PI,kr);
  Kcomp_th = ellK(0.5*PI,kth);
   
/*---Computing Incomplete elliptical integrals at r(lambda=0) = rperi---*/
  yr = 0.0;
  yth = 0.0;
  Fin_r = ellK(yr,kr);
  Fin_th = ellK(yth,kth);
    
/*---Computing radial and polar times at r(0)=rperi---*/
  epsilon0 = SQR(GV_Spin)*(1.0 - SQR(GV_E))/SQR(GV_Lz);
    
  lambda1_r = 2.0*Fin_r/(sqrt((1.0-SQR(GV_E))*(GV_R_apo -GV_R_3)*(GV_R_peri-GV_R_4)));
  lambda1_th = Fin_th/(GV_Lz*sqrt(epsilon0*GV_Z_plus));
   
/*---Computing the analytical radial solution---*/
/**---Computing Jacobi's elliptic function----**/
  M_r1 = lambda_r/2 - lambda1_r;
  M_r2 = lambda_r - lambda1_r;
    
  if( (dMt_r >= 0.0)||(dMt_r <= M_r1) )
  {printf("hola r1\n");
    ur_M = dMt_r + lambda1_r;
    ur_M /= lambda_r;
    ur_M *= 2.0*Kcomp_r;
  }
  if( (dMt_r >= M_r1 )||(dMt_r <= M_r2) )
  {printf("hola r2\n");
    ur_M = -dMt_r + lambda_r - lambda1_r;
    ur_M /= lambda_r;
    ur_M *= 2.0*Kcomp_r;

  }
  if( (dMt_r >= M_r2 )||(dMt_r <= lambda_r) )
  {printf("hola r3\n");
    ur_M = dMt_r - lambda_r + lambda1_r;
    ur_M /=lambda_r;
    ur_M *= 2.0*Kcomp_r;
  }
  else exit(0);
    
  //  printf("dMt= %4.6e, dMt_r = %4.6e, M_r1  = %4.6e, M_r2 = %4.6e, N = %d  \n\n",dMt, dMt_r, M_r1, M_r2,omega_r, floor(omega_r*dMt/(2.0*PI)) );

   JacobianEllipticFunctions(ur_M,kr, &sn_ur, &cn, &dn);
      
   r_M = (GV_R_3*(GV_R_apo - GV_R_peri)*SQR(sn_ur) - GV_R_peri*(GV_R_apo-GV_R_3))/((GV_R_apo-GV_R_peri)*SQR(sn_ur) - (GV_R_apo-GV_R_3));

   GV_Rres[l_time] = r_M;
     
/*---Computing the analytical polar solution---*/
/**---Computing Jacobi's elliptic function----**/
  M_th1 = lambda_th/4.0 - lambda1_th;
  M_th2 = 3.0*lambda_th/4.0  - lambda1_th;
    
  if( (dMt_th >= 0.0)||(dMt_th <= M_th1) )
  {
      printf("hola th1\n");
      uth_M = dMt_th + lambda1_th;
      uth_M /= lambda_th;
      uth_M *= 4.0*Kcomp_th;
  }
  if( (dMt_th >= M_th1 )||(dMt_th <= M_th2) )
  {
      uth_M = -dMt_th + 0.5*lambda_th - lambda1_th;
      uth_M /= lambda_th;
      uth_M *= 4.0*Kcomp_th;
      printf("hola th2\n");
  }
  if( (dMt_th >= M_th2 )||(dMt_th <= lambda_th) )
  {
      printf("hola th3\n");
      uth_M = dMt_th -lambda_th + lambda1_th;
      uth_M /= lambda_th;
      uth_M *= 4.0*Kcomp_th;
  }
  else exit(0);
   // printf(" dMt_th = %4.6e,  M_th1  = %4.6e, M_th2 = %4.6e, omega_th = %4.6e\n",dMt_th, M_th1, M_th2, omega_th );
   
 
  JacobianEllipticFunctions(uth_M, kth, &sn_uth, &cn, &dn);
   
  theta_M = sqrt(GV_Z_minus)*sn_uth;
  theta_M = acos(theta_M);
    
  GV_Thres[l_time] = theta_M;
    
/*----Creating directory to store analytical solutions----*/
    strcpy(GV_Fn_analytical_solutions, GV_Dn_run);
    strcat(GV_Fn_analytical_solutions, "/Analytical_Solutions_R");
   
  data = fopen(GV_Fn_analytical_solutions,"a");
   fprintf(data," %4.6e	 %4.6e	%4.6e   %4.6e	%4.6e \n", dMt, dMt_r, dMt_th, GV_Rres[l_time], GV_Thres[l_time]);
 fclose(data);
   
/*----- This is the end -----*/
    return 0;
}








