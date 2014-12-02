/*==============================================================================
                            GLOBAL VARIABLES
 
        NOTE:  The Mass of the Black Hole has been set equal to one.
               Global quantities convention si that they start by GV_
 ------------------------------------------------------------------------------
                                            Created by Priscilla on 19/12/2012
 ------------------------------------------------------------------------------
  Last Update : 30.12.12
 =============================================================================*/


/*-----  MATHEMATICAL LIBRARIES  -----*/
/**----- INTEL math library -----**/
#ifdef INTEL_MATH
#include <mathimf.h>
#endif

/**----- GNU math library -----**/
#ifdef GNU_MATH
#include <math.h> 
#endif  


/*----- Algebraic quantities ----*/
  double **Kronecker;         // Spatial Kronecker delta (it has four components in order to use indices from 1 to 3 instead that from 0 to 2)
  double Levi_Civita[4][4][4];// Spatial Levi-Civita antisymmetric symbol (it has four components in order to use indices from 1 to 3 instead that from 0 to 2)

/*----- EMRI PARAMETERS -----*/
  double EMRI_Parameter[18];            // Parameters that characterize the EMRI
  double PSI_o;                         // Initial values of Psi_p_o.
  double CHI_o;                         // Initial values of Chi_p_o.
  double PHI_o;                         // Initial values of PHI_o.
  double R_p;                           // Radial coordinate .
  double COSTh_p;                       // Polar coordinate.
  double TAU, TIME;                     // Proper & coordinate times
  double Qr, Qth;                       // Angle variables
  double MAG_SF;                        // Sefl-force magnitude
  double Lambda,Lambdar,Lambdath;

/*----- WAVEFORM PARAMETERS -----*/
 double **E_p, **E_x;                  // Waveform polarization tensors
 double dRdt, dCOSTHdt;
 double GAMMA[4][4][4];                // Christoffel symbols of the Kerr metric associated with Boyer-Lindquist

/*----- COMPUTATIONAL PARAMETERS -----*/
  double Real_Parameter[18];             // (Real) Computational Parameters, diverse types
  long unsigned Integer_Parameter[18];   // (Integer) Computational Parameters, diverse types
  
/*----- Strings for Directory Names -----*/
  char DirRun[100];                      // Directory to store the data of a  Run
  char FN_SCO_evolution[100];
  char FN_rva_particle[100];		     // FIle to store the  position, velocity and acceleration 
  char FN_frequencies[100];              // File to store the fundamental frequencies associated with a geodesic orbit
  char FN_Res_frequencies[100];          // File to store the fundamental frequencies on resonance
  char FN_geo_evo[150];                  // File to store the geodesic evolution
  char FN_parameters[150];               // File to store run parameters
  char FN_evolution_AA[150];             // File to store the osculating evolution
  char FN_evolution_RR[150];             // File to store the radiative evolution
  char FN_fluxes[150];                   // File to store the fluxes of E, lz and Q
  char FN_waveforms[150];                // File to store the waveforms

/*----- Strings for File Names -----*/
  char FN_LogFile[100];                  // Log File where the main activity of the Code is written

/*----- Flag that control outputs  -----*/
  char FLAG_RR;                          // Flag that controls whether we evolve the orbits using Radiative effects
  char FLAG_OS;                          //  " use action angle variable formalism
  char FLAG_resonances;                  // Flag to Choose whether we consider resonances
  char FLAG_incConsPiece;                // parameter to control the inclusion of the conservative piece of the Self-force

/*---- Self-Force Parameters----*/
  double EPSILON;                         // Force magnitude
  double PNFcomp[4];                      // Post-Newtonian Self-Force components t =[0], r =[1], theta = [2], phi = [3]
  double STORE_rate;                      // Frequency for data saving
