/*==============================================================================
 
        This routine initializes the variables used in the current run
 
 -------------------------------------------------------------------------------
                                             Created by Priscilla on 20/12/2012
 -------------------------------------------------------------------------------
 Last Update: 25.09.13
 =============================================================================*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "physical_quantities.h"
#include "global_quantities_osc.h"

#include "macros_2.h" 

int checks_errors(int);
int compute_Keplerian_period(double, double, double *, long unsigned *, double);
int pez_to_ELCQ(char, double, double, double,double, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int turning_points(double, double, double, double, double, double, double *, double *, double *, double *, double *);
int compute_frequencies(double *, double *, double *, double *, double *,double *, double *, double *, double *);
int scalar_product(double *, double *, double *);
int vector_product(double *, double *, double *);
int vector_norm(double *, double *);

double *Realvector(const long , const long );
double **Realmatrix(const long , const long , const long , const long );

void free_Realvector(double *, const long, const long);

int initialization(void)
{
  double *sbs, *source, *spinvec;	// solar base system vector, source vector, spin vector
  double norm;
  double *p_src, *q_src;		// Normal vectors (in the Solar Barycentric System) that form, together with n, a spatial orthonormal triad with respect to the Euclidean scalar product defined via the Kronecker delta (independently of the Cartesian coordinates choosen)
  int i,j,k;

/*----- Setting Global System Parameters for this Run -----*/
  MASS_in_MSUN = EMRI_Parameter[1];          // Black Hole Mass (in units of the Mass of the Sun)
  SPIN         = EMRI_Parameter[2];          // Black Hole Spin (in units of the Black Hole Mass)
  MASS_SCO     = EMRI_Parameter[3];          // Mass of the Particle orbiting the Black Hole (in units of the Black Hole Mass)
  ECCENTRICITY = EMRI_Parameter[4];          // Initial Eccentricity
  P_p          = EMRI_Parameter[5];          // Initial Semilatus Rectum
  // See conventions in PRD 86, 044010 (2012)
  THETA_inc    = EMRI_Parameter[6]*PI/180.0; // Initial Inclination of the Orbit in radians
  THETA_Source = EMRI_Parameter[7];          // EMRI Polar Angle with respect to the Ecliptic Reference Frame.
  PHI_Source   = EMRI_Parameter[8];          // EMRI Azimuthal Angle with respect to the Ecliptic Reference Frame.
  THETA_Kerr   = EMRI_Parameter[9];          // Spin Vector Polar Angle with respect to the Ecliptic Reference Frame.
  PHI_Kerr     = EMRI_Parameter[10];         // Spin Vector Azimuthal Angle with respect to the Ecliptic Reference Frame.
  DL           = EMRI_Parameter[11];	     // Luminosity Distance divided by the SCO mass.
  PSI_o        = EMRI_Parameter[12];         // Boyer-Lindquist Angle psi [associated with r_BL] at t = t_o: psi_o.
  CHI_o        = EMRI_Parameter[13];         // Boyer-Lindquist Angle chi [associated with theta_BL] at t = t_o: chi_o.
  PHI_o        = EMRI_Parameter[14];         // Boyer-Lindquist Angle varphi [identical to varphi_BL] at t = t_o: varphi_o.

/*----- Setting Global computational Parameters for this Run -----*/
  EVAL_arg = 0;                             // Counter to track whether the elliptical integrals have been passed invalid arguments
  DT_MBH = 10.0*DT_INI;                     // ODE Global Time Step
  N_periods = 1.;                           // Number of "periods" between evolving a geodesic

/*---- Setting Flags for this Run -----*/
 if ( SPIN < 1.0e-14 )
  {
    FLAG_Schwarzschild = ‘y’;		// Flag that controls whether the Black Holes is described by the Schwarzschild metric tensor.  VALUES: y (yes), n (no)
    SPIN = 0.0;
  }
  else FLAG_Schwarzschild = 'n';

  if (FLAG_Schwarzschild == 'y')
  {
    Q_CONSTANT= 0.0;
    
    FLAG_equatorial_orbit = 'y';
    THETA_inc = 0.0;
  }


  if ( ECCENTRICITY < 1.0e-14 )
  {
    FLAG_CircOrb = 'y';                // Flag that controls whether the orbit is circular
    ECCENTRICITY = 0.0;
  }
  else FLAG_CircOrb = 'n';
  
  if ( THETA_inc < 1.0e-14 )
  {
    FLAG_equatorial_orbit = 'y';       // Flag that controls whether the orbit is equatorial or not
    THETA_inc = 0.0;
  }
  else FLAG_equatorial_orbit = 'n';
  
  if ( THETA_inc < 0.0 )             // Flag that controls whether the orbit is prograde or retrograde
    FLAG_orbit_spin = 'r';
  else 
    FLAG_orbit_spin = 'p';           
/**----- Inclination Angle must be in the interval [-PI/2,PI/2] -----**/
    if ( (THETA_inc < -0.5*PI) || (THETA_inc > 0.5*PI) )  checks_errors(3);

//----Orbital Phases initialization
  if(FLAG_RR=='y')
  {
    if(FLAG_equatorial_orbit == 'y' )
    {
      PSI_OLD = PSI_o;
      CHI_OLD = 0.5*PI;
      PHI_OLD = PHI_o;
    }
    else
    {
      PSI_OLD = PSI_o;
      CHI_OLD = CHI_o;
      PHI_OLD = PHI_o;
    }
  }
  else 
  {
    PSI_OLD =  0.0;
    CHI_OLD =  0.0;
    PHI_OLD =  0.0;
    
    TAU = 0.0;   // Proper time
    TIME = 0.0;  // coordinate time
    Qr  = 0.00;
    Qth = PI/2.;
    ODEcount = 0;
    Lambda = 0.0;
    
    //  ENERGY_dot = 0.0;
    //LZ_dot = 0.0;
    //Q_CONSTANT_dot = 0.0;
  }
  
/*----- Global Black Hole Parameters [M^2,S^2,R_+,R_-,sqrt(M^2-a^2)] -----*/
  MASS_ratio = MASS_SCO;                    // mu = m/M
  MASS_ratio2 = SQR(MASS_ratio);            // mu^2
  SPIN2  = SQR(SPIN);                       // (a/M)^2
  R_plus = 1.0 + sqrt(1.0 - SPIN2);          // r+ = M + sqrt(M^2 - a^2)
  R_minus= 1.0 - sqrt(1.0 - SPIN2);          // r- = M - sqrt(M^2 - a^2)
 
/*-----Computing Inital Constants of Motion -----*/
  pez_to_ELCQ(FLAG_CircOrb, P_p, ECCENTRICITY, THETA_inc,SPIN,
              &INCLINATION, &R_peri, &R_apo, &THETA_min, &Z_minus,
              &ENERGY, &LZ, &C_CONSTANT, &Q_CONSTANT);
    
/*-----Computing initial turning points of the Boyer-Linquist radial coordinate 'r' and 'theta' -----*/
  turning_points(P_p, ECCENTRICITY, ENERGY, LZ,Q_CONSTANT,Z_minus, &R_3, &R_4, &P3_p, &P4_p, &Z_plus);

/*----Computing Initial fundamental frequencies (associated with Mino time) and initial arguments of the Elliptic functions ---*/
  compute_frequencies(&OMEGA_r, &OMEGA_th, &OMEGA_t, &Kcom_r, &Kcom_th, &Ecom_r, &Ecom_th, &KR, &KTH);

/*----- Computing Keplerian period and number of time steps inside it (used here to apply RR) -----*/
  if(FLAG_RR=='y')
   compute_Keplerian_period(P_p,ECCENTRICITY,&KEPLER_period,&StepKT,DT_MBH);
  if ( (FLAG_RR == 'y') && ((StepKT*N_periods) > NMAX) ) checks_errors(23);

//=======  WAVEFORMS VARIABLES =====//
  
/*----- Setting Time Independent Quantities for GW observations -----*/
  
  DL = DL*(1.0e9*PARSEC)/BHMASS_to_m; // Distance from the Solar System Barycenter to the Source/EMRI (D_L) in MBH mass

/**----- Unit Vector that points from the detector to the GW source written in the SBS reference frame [NOTE: In general this could be a time dependent vector] -----**/
  sbs = Realvector(1, 4);
  sbs[1] = sin(THETA_Source)*cos(PHI_Source);
  sbs[2] = sin(THETA_Source)*sin(PHI_Source);
  sbs[3] = cos(THETA_Source);
    
/**----- Unit Vector that points from LISA to the EMRI written in the Source reference frame [NOTE: This is an approximation.  In general this is a time dependent vector] -----**/
  source = Realvector(1, 4);
  source[1] = cos(THETA_Kerr)*cos(PHI_Kerr)*sbs[1] + cos(THETA_Kerr)*sin(PHI_Kerr)*sbs[2] - sin(THETA_Kerr)*sbs[3];
  source[2] = -sin(PHI_Kerr)*sbs[1] + cos(PHI_Kerr)*sbs[2];
  source[3] = sin(THETA_Kerr)*cos(PHI_Kerr)*sbs[1] + sin(THETA_Kerr)*sin(PHI_Kerr)*sbs[2] + cos(THETA_Kerr)*sbs[3];
    
/**----- Unit Vector that Points in the MBH Spin written in the Source reference frame -----**/
  spinvec = Realvector(1, 4);
  spinvec[1] = 0.0;
  spinvec[2] = 0.0;
  spinvec[3] = 1.0;
    
/*----- Computing the unit vectors p and q -----*/
/*----- n x s -----*/
 p_src = Realvector(1, 4);
 vector_product(p_src, source, spinvec);
    
/*----- | n x s | -----*/
 vector_norm(&norm, p_src);
    
/*----- p = ( n x s ) / | n x s | -----*/
   for (i=1;i<=3;i++)
    p_src[i] /= norm;
    
/*----- q = p x n -----*/
 q_src = Realvector(1, 4);
 vector_product(q_src, p_src, source);
  

/*----- Computing Polaritzation tensors e_p, e_x [See conventions in PRD 86, 044010 (2012)] -----*/
  E_p = Realmatrix(1, 4, 1, 4);
  E_x = Realmatrix(1, 4, 1, 4);
    
    for (i=1;i<=3;i++)
    {
        for (j=1;j<=3;j++)
        {
            E_p[i][j] = (p_src[i])*(p_src[j]) - (q_src[i])*(q_src[j]);
            E_x[i][j] = (p_src[i])*(q_src[j]) + (q_src[i])*(p_src[j]);
          
        }
    }
  
  free_Realvector(p_src, 1, 4);
  free_Realvector(q_src, 1, 4);
  
  free_Realvector(source, 1, 4);
  free_Realvector(spinvec,1, 4);
  free_Realvector(sbs,1, 4);
  
/*----- This is the end -----*/
  return 0;
}

