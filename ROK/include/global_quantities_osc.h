//
//  global_quantities_osc.h
//  X_ROK_V1.6
//
//  Created by Priscilla on 26/04/2013.
//
//

/*-----  MATHEMATICAL LIBRARIES  -----*/
/**----- INTEL math library -----**/
#ifdef INTEL_MATH
#include <mathimf.h>
#endif

/**----- GNU math library -----**/
#ifdef GNU_MATH
#include <math.h>
#endif

/*----- Resonace Study -----*/
double OMEGA_r, OMEGA_th;          // Radial and Polar frequencies in "Mino" time
double OMEGA_phi, OMEGA_t;         // Azimuthal and temporal frequencies in "Mino" time
double OMEGA_r_res, OMEGA_th_res;  // Radial and Polar frequencies  "in Mino time", on resonance

//double Ec_o_search, P_o_search;    // Initial values of the eccentricity and semilatus rectum in resonance search
//double EC_f_search, P_f_search;    // Final values of the eccentricity and semilatus rectum in resonance search
//double DEC, DP;                    // Eccentricity and semi-latus rectum step size
double Res_accuracy;               // Resonant condition (k*w_th - n*w_r) accuracy

double Kcom_r, Kcom_th;            // Complete elliptic integral of the firs kind
double Ecom_r, Ecom_th;            // Complete elliptic integral of the second kind
double KR, KTH;                    // Arguments of the complete elliptic integrals
                                   
int Resorder[3];                   // Resonance order (k*w_r - n*w_th): Resorder[0] = k, Resorder[1] = n
int RES_COUNT;                     // Counter to track resonances
int NVAR;                          // Number of variables evolved in the osculating method
int ODEcount;