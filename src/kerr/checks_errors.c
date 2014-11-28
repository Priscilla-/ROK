/*=====================================================================================
 
   This Routine checks consistency and compatibility of different parts of the Code.
 
 ---------------------------------------------------------------------------------------
   Last Update : 28.12.12   by PCM
 ======================================================================================*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include "global_quantities_kerr.h"

int checks_errors(int type)
{

/*----- Switch that selects the particular check of a given call -----*/
  switch(type)
  {
    case 1 :
             /*----- Shows how to Run the Compiled ROK Code -----*/
	         printf("  |--------------------------------------------------------------------------------\n  | \n");
             printf("  | Usage: ./exe/ROK Run_Number\n  | \n");
             printf("  | Format of Run_Number: XXXX (Four Digits!)\n  | \n");
             printf("  |--------------------------------------------------------------------------------\n  |\n");
             exit(1);
	        
    case 2 :
             /*-----  -----*/
	         printf("  |--------------------------------------------------------------------------------\n  | \n");
	         printf("  |--------------------------------------------------------------------------------\n\n");
             exit(1);
	        
    case 3 :
             /*----- Inclination Angle in Parameter File out of range -----*/
	         printf("  |--------------------------------------------------------------------------------\n  |\n");
			 printf("  | The inclination angle specified in the Parameter File is:\n  |\n");
			 printf("  |            Theta_inc = %4.8f ,\n  |\n",THETA_inc);
			 printf("  | therefore it is out of its range: [-PI/2,PI/2]. Please, change this and try again.\n  |\n");
	         printf("  |--------------------------------------------------------------------------------\n");
             exit(1);
	        
    case 4 :
             /*----- Something Wrong has happened in the Routine 'pez_to_ELCQ.c' -----*/
	         printf("  |--------------------------------------------------------------------------------\n  | \n");
             printf("  | There must some inconsistency in the computation of the constants of motion from\n");
			 printf("  | the orbital initial data (p,e,theta_inc).  This has been detected by the diagnostics\n");
			 printf("  | of the routine 'pez_to_ELCQ.c'.  Please, check this routine to locate the problem.\n");
             printf("  |--------------------------------------------------------------------------------\n\n");
             exit(1);
	        
    case 5 :
			 /*-----  -----*/
	         printf("  |--------------------------------------------------------------------------------\n  | \n");
             printf("  |--------------------------------------------------------------------------------\n\n");
             exit(1);
	        
    case 6 :
             /*----- ERROR IN MEMORY ALLOCATION -----*/
	         printf("  |--------------------------------------------------------------------------------\n  | \n");
             printf("  | ERROR allocating Memory.  The Code could not allocate the memory\n");
             printf("  | that you asked for and has to exit.  Sorry about that! Goodbye!\n  |\n");
	         printf("  |--------------------------------------------------------------------------------\n\n");
             exit(1);

    case 7 :
             /*-----  -----*/
	         printf("  |--------------------------------------------------------------------------------\n  | \n");
             printf("  |--------------------------------------------------------------------------------\n\n");
             exit(1);
	        
    case 8 :
          /*-----  -----*/
          printf("  |--------------------------------------------------------------------------------\n  |\n");
          printf("  |--------------------------------------------------------------------------------\n\n");
          exit(1);
	case 9:
	         /*-----  -----*/
			 printf("  |--------------------------------------------------------------------------------\n  |\n");
             printf("  |--------------------------------------------------------------------------------\n\n");
             exit(1);
			 	        
	case 10:
	         /*----- Feature not implemented -----*/
			 printf("  |--------------------------------------------------------------------------------\n  |\n");
			 printf("  | One of the features of the Code that you have selected has not been implemented yet.\n");
			 printf("  | Sorry about that.  Goodbye\n  |\n");
             printf("  |--------------------------------------------------------------------------------\n\n");
             exit(1);
			 	        
    case 11 :
             /*-----  -----*/
	         printf("  |--------------------------------------------------------------------------------\n  | \n");
             printf("  |--------------------------------------------------------------------------------\n");
             exit(1);
	        
	case 12:
	         /*-----  -----*/
			 printf("  |--------------------------------------------------------------------------------\n  |\n");
             printf("  |--------------------------------------------------------------------------------\n\n");
             exit(1);
			 	        
    case 13 :
             /*----- Wrong Flag for the Spin of the Orbit -----*/
	         printf("  |--------------------------------------------------------------------------------\n  |\n");
			 printf("  | The flag for the 'Spin' of the Orbit has a wrong value:\n  |\n");
			 printf("  |            FLAG_orbit_spin = %c ,\n  |\n",FLAG_orbit_spin);
			 printf("  | it can only take the values 'p' (prograde orbits) and 'r' (retrograde orbits).\n  |\n");
			 printf("  | Please, change this and try again.\n  |\n");
	         printf("  |--------------------------------------------------------------------------------\n");
             exit(1);
	        
    case 14 :
             /*----- Incompatibility of the Inclination angle with the Flag for the Spin of the Orbit -----*/
	         printf("  |--------------------------------------------------------------------------------\n  |\n");
			 printf("  | The flag for the 'Spin' of the Orbit and the Inclination angle (Theta_inc) have values\n");
			 printf("  | within their ranges but they are incompatible.  That is, either Theta_inc is less than\n");
			 printf("  | PI/2 and the orbit is retrograde, or Theta_inc is greater than PI/2 and the orbit is\n");
			 printf("  | prograde.  Please, change this in the parameter file 'parameters_evolution'.\n  |\n");
	         printf("  |--------------------------------------------------------------------------------\n");
             exit(1);

	case 15:
	         /*----- Incompatibility between the flag 'FLAG_circular_orbit" and the value of the Eccentricity -----*/
			 printf("  |--------------------------------------------------------------------------------\n  |\n");
			 printf("  | The flag 'FLAG_circular_orbit' was set equals to 'y'.  However, the value of the\n");
			 printf("  | Eccentricity is:\n");
			 printf("  |                      Eccentricity = %2.8e\n  |\n", ECCENTRICITY);
			 printf("  | which for the standards of this code is by no means zero as it should be.  Please,\n");
			 printf("  | either change the flag or the value of the Eccentricity, and then run the code again.\n");
			 printf("  | Goodbye!\n  |\n");
             printf("  |--------------------------------------------------------------------------------\n");
             exit(1);

    case 16 :
	         /*----- Incompatibility between the flag 'FLAG_equatorial_orbit" and the value of the inclination angle Theta_inc -----*/
	         printf("  |--------------------------------------------------------------------------------\n  |\n");
			 printf("  | The flag 'FLAG_equatorial_orbit' was set equals to 'y'.  However, the value of the\n");
			 printf("  | inclination angle Theta_inc is:\n");
			 printf("  |                      Theta_inc = %2.8e\n  |\n", THETA_inc);
			 printf("  | which for the standards of this code is by no means zero as it should be.  Please,\n");
			 printf("  | either change the flag or the value of Theta_inc, and then run the code again.\n");
			 printf("  | Goodbye!\n  |\n");
             printf("  |--------------------------------------------------------------------------------\n");
             exit(1);
	        
    case 17 :
             /*----- Bad Selection of the Flag that controls whether the Spacetime is Schwarzschild or that of a Spinning Black Hole -----*/
	         printf("  |--------------------------------------------------------------------------------\n  |\n");
             printf("  | The Flag 'FLAG_Schwarzschild' has been set up to the value:\n  |\n");
			 printf("  |       FLAG_Schwarzschild = '%c'\n  |\n",FLAG_Schwarzschild);
			 printf("  | However, the only valid values are either 'y' or 'n'. Please, change the parameter\n");
			 printf("  | file to fix this.  Goodbye!\n  |\n");
             printf("  |--------------------------------------------------------------------------------\n");
             exit(1);
			 
    case 18 :
             /*----- Schwarzschild Black Hole should be chosen with Equatorials Orbits -----*/
	         printf("  |--------------------------------------------------------------------------------\n  |\n");
             printf("  | The Flag 'FLAG_Schwarzschild' has been set up.  The most convenient thing for such\n");
			 printf("  | a type of run is to set up also the Flag 'FLAG_equatorial_orbit'. Please, change\n");
			 printf("  | the parameter file to fix this.  Goodbye!\n  |\n");
             printf("  |--------------------------------------------------------------------------------\n");
             exit(1);
			 
    case 19 :
			 /*----- There are not anymore 4 real roots in the Radial Motion -----*/
			 printf("  |--------------------------------------------------------------------------------\n  |\n");
			 printf("  |  The Code has detected that there are no anymore 4 REAL Roots in the Radial\n");
			 printf("  |  Motion [Roots of the Polynomial R(r)].\n");
			 printf("  |  Therefore, the SCO has plunged into the MBH!!!\n  |\n");
             printf("  |--------------------------------------------------------------------------------\n");
			 exit(1);

    case 20 :
			 /*----- Plunge according to Criteria establish in routine 'check_if_orbit_is_plunge.c' -----*/
			 printf("  |--------------------------------------------------------------------------------\n  |\n");
			 printf("  |  According to some criteria establish in this code, the orbit, after applying\n");
			 printf("  |  Radiation Reaction Effects has become a plunge Orbit.\n");
			 printf("  |  Therefore, the SCO has plunged into the MBH!!!\n  |\n");
             printf("  |--------------------------------------------------------------------------------\n");
			 exit(1);

    case 21 :
			 /*-----  -----*/
			 printf("  |--------------------------------------------------------------------------------\n  |\n");
             printf("  |--------------------------------------------------------------------------------\n");
			 exit(1);
			 
    case 22 :
			 /*-----  -----*/
			 printf("  |--------------------------------------------------------------------------------\n  |\n");
             printf("  |--------------------------------------------------------------------------------\n");
			 exit(1);
			 
    case 23 :
			 /*----- Incompatibility in the Number of Time Steps -----*/
			 printf("  |--------------------------------------------------------------------------------\n  |\n");
			 printf("  |  You want to perform an evolution that takes into account Radiation Reaction\n");
			 printf("  |  effects, but the number of Total Time Steps is:\n  |\n");
			 printf("  |              Total Number of Time Steps = %lu\n  |\n",NMAX);
			 printf("  |  which is smaller than the number of times steps that the code will evolve\n");
			 printf("  |  between the application of Radiation Reaction effects, that is:\n  |\n");
			 printf("  |   # Time Steps without RR = (# Kepler Periods) x (# Time Steps per Kepler Period)\n");
			 printf("  |                           =        %lu  x  %lu  =  %lu .\n  |\n",N_periods,StepKT,N_periods*StepKT);
			 printf("  |  Then, Radiation Reaction will never be used in this Run. Please, change\n");
			 printf("  |  the parameter file to fix this.  Goodbye!\n  |\n");
             printf("  |--------------------------------------------------------------------------------\n");
			 exit(1);
			 
    default :
             printf("  | The \"type\" number must be wrong.  Check this.  Now the code has\n");
             printf("  | to stop. Goodbye!\n  | \n");
	         printf("  |--------------------------------------------------------------------------------\n\n");
             exit(1);
             break;
  }


/*----- This Routine ends Here -----*/
  return 0;

}

