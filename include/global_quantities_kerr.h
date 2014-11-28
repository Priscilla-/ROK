/*==============================================================================
                            GLOBAL VARIABLES
 
   NOTE:  The Mass of the Black Hole has been set equal to one.  
          Global quantities convention si that they start by GV_
 ------------------------------------------------------------------------------
                                            Created by Priscilla on 19/12/2012
 ------------------------------------------------------------------------------
 Last Update : 11.12.12                                                                                       
==============================================================================*/


/*-----  MATHEMATICAL LIBRARIES  -----*/
/**----- INTEL math library -----**/
#ifdef INTEL_MATH 
#include <mathimf.h>                             
#endif 

/**----- GNU math library -----**/
#ifdef GNU_MATH
#include <math.h>                                       
#endif


/*-----  UNIVERSAL VARIABLES  -----*/
#define PI 3.1415926535897932384626433832795029
  
/*----- BLACK HOLE QUANTITIES -----*/
  double MASS_in_MSUN;                               // Mass of the Black Hole in Solar Masses
  double SPIN, SPIN2;                             // Spin of the Black Hole (in units of the Black Hole Mass)
  double R_plus, R_minus;                         // Black Hole horizon location in Boyer-Lindquist coordinates.  R_plus is the true horizon, whereas R_minus is always inside R_plus (R_minus <= R_plus)

/*-----	PARTICLE PROPERTIES -----*/
  double MASS_SCO;                                   // Mass of the Particle orbiting the Black Hole  (in units of the Black Hole Mass)
  double MASS_ratio;                                 // Mass Ratio: Mass_ratio = Mass_SCO 


/*----- OTHER EMRI PARAMETERS -----*/
  double THETA_Source;                               // EMRI Polar Angle with respect to the Ecliptic Reference Frame.
  double PHI_Source;                                 // EMRI Azimuthal Angle with respect to the Ecliptic Reference Frame.
  double THETA_Kerr;                                 // Spin Vector Polar Angle with respect to the Ecliptic Reference Frame.
  double PHI_Kerr;                                   // Spin Vector Azimuthal Angle with respect to the Ecliptic Reference Frame.
     
  double ECCENTRICITY;                               // Eccentricity of the Orbit
  double P_p;                                        // Semi-latus rectum of the Orbit  
  double INCLINATION;                                // Inclination of the Orbit
  double THETA_inc;                                  // Another variable for the Inclination of the Orbit
  double DL;                                         // Luminosity Distance to the EMRI from the Ecliptic Baricentric Reference Frame.

  double ENERGY;                                      // Energy of the Orbit  
  double LZ;                                          // Angular Momentum of the Orbit  
  double Q_CONSTANT,C_CONSTANT;                       // Carter constant (two versions of it).  They are related by: Q = C + ( a E - Lz )^2.  Have a look at the file 'turning_points.c' to see how they appear in the geodesic equations


  double R_apo,R_peri,R_3,R_4;              // Turning points of the 'radial' Boyer-Lindquist coordinate 'r'
  double THETA_min;                         // Minimum value of Theta
  double Z_minus,Z_plus;                    // Turning points of the Boyer-Lindquist angular coordinate 'theta'

/*----- Useful Auxiliary Quantities -----*/  
  double MASS_ratio2;                                // Auxiliary quantities to make more efficient the computation of the ODE RHSs
  double P3_p,P4_p;                                  //

/*----- Some Variables used for the Radiation Reaction Computations -----*/
  long unsigned StepKT;                             // Number of Time Steps in a single (initial) Kepler Period
  long unsigned N_periods;                          // Number of Kepler Periods to be evolved without introducing RR effects
  long unsigned N_one_Radial_period;                // Number of Time Steps in a single (initial) Radial Period

  double KEPLER_period;                             // Keplerian Period associated with the orbital parameters in GR
  double ENERGY_dot;                                // Time Evolution of the Energy of the Orbit
  double LZ_dot;                                    // Time Evolution of the Angular Momentum in the Spin direction of the Orbit
  double Q_CONSTANT_dot;                            // Time Evolution of the Carter Constant Q_carter of the Orbit


/*----- ORBITAL QUANTITIES USED IN THE EVOLUTION ALGORITHM -----*/
  double PSI_OLD, PSI_NEW;                    // ODE Coordinate Psi associated with the Boyer-Lindquist Radial Coordinate R_p in order to avoid turning points in the motion
  double CHI_OLD, CHI_NEW;                    // ODE Coordinate Chi associated with the Boyer-Lindquist Polar Coordinate Theta_p in order to avoid the turning points in the motion
  double PHI_OLD, PHI_NEW;                    // Boyer-Lindquist Azimutal Coordinate Phi_p
  

 double R_p_o, COSTH_p_o;                       // We only need to add this two with respect to the ODE coordinate system

 double R_p_oLD , R_p_NEW;
 double COSTH_p_oLD , C_CONSTANT;
 

/*----- ODE Integrator Parameters -----*/
  double TIME_yrs;                          // Evolution Time [yrs]
  double DTsec;                             // Sampling Time/Data Saving Time Interval [sec]
  double DT_MBH;                            // ODE Global Time Step
  double DT_INI, DT_MIN, ACCURACY;          // Time Steps (initial and minimum) in Black Hole Mass and accuracy level for ODE integrator 
  
  long unsigned NMAX, Nfreq;                // Number of Time steps and Time Frequencies of Data Saving
  long unsigned NSAMPLING;                  // Total Number of Sampling Times in the Time Series that the Code uses


/*----- Flags and Control Parameters -----*/
  char FLAG_circular_orbit;                          // Flag that controls whether the orbit is circular
  char FLAG_orbit_spin;                              // Flag that controls whether the orbit is prograde or retrograde.  VALUES: p (prograde), r (retrograde) 
  char FLAG_Schwarzschild;                           // Flag that controls whether the Black Holes is described by the Schwarzschild metric tensor.  VALUES: y (yes), n (no)
  char FLAG_equatorial_orbit;                        // Flag that indicates that the orbits is essentially equatorial, so we can assume that Theta_p = Pi/2 and DTheta_p_DT = 0.  VALUES: y (yes), n (no)

/*----- Strings for File Names -----*/	  
  char FN_LogFile[100];                             // Log File where the main activity of the Code is written
  char FN_analytical_solutions[150];                // File to store the analytical solutions

/*---- Ellipitical ----*/
  int EVAL_arg;                                    // Counter to determine whether the elliptical integrals have been passed invalid arguments








