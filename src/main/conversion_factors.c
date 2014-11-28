/*==============================================================================
 
     This file computes some conversion factors that are useful to convert the 
     results into meaningful and/or useful quantities.
 -------------------------------------------------------------------------------
                                             Created by Priscilla on 19/12/2012
 -------------------------------------------------------------------------------
 Last Update : 6.07.13
 ==============================================================================*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "global_quantities_main.h"
#include "physical_quantities.h"
#include "macros_2.h"


int conversion_factors(void)
{

/*----- Units, Physical Constants, and Conversion Factors -----*/
  CSPEED3 = SQR(LIGHT_SPEED)*LIGHT_SPEED;       // c^3

  BHMASS_to_s = (MASS_SUN*EMRI_Parameter[1])*(GN/CSPEED3);           // Factor to convert Times in MBH mass to Times in Seconds
  
  BHMASS_to_m = (GN*EMRI_Parameter[1])*(MASS_SUN/SQR(LIGHT_SPEED));  // Factor to convert Distances in MBH mass to Distances in Meters
    
    
  //GC_Conversion_from_meters_to_parsecs = 1.0/PARSEC;         // Factor to convert lengths in meters to parsecs
  // GC_Conversion_from_seconds_to_years = 1.0/SECSYR;          // Factor to convert times in secons to sideral years

    
/*----- this is the end -----*/
  return 0;
}
