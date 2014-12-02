/*=======================================================================================
                                                                                       
   PHYSICAL QUANTITIES:  As the Global quantities, they are visible from everywhere in 
   the code (the convention is that they should always start by a capital letter.
 ---------------------------------------------------------------------------------------
  Last Update : 23.08.12                                                               
==========================================================================================*/


/*-----  MATHEMATICAL CONSTANTS -----*/
#define PI 3.1415926535897932384626433832795029

/*----- NEWTON'S CONSTANT [m^3 Kg^-1 s^-2] -----*/
// The gravitational constant, denoted G, is an empirical physical constant involved in 
// the calculation of the gravitational attraction between objects with mass. It appears 
// in Newton's law of universal gravitation and in Einstein's theory of general relativity. 
// It is also known as the universal gravitational constant, Newton's constant, and 
// colloquially Big G.

#define GN 6.67428e-11


/*----- SIDEREAL YEAR [s] -----*/
// The sidereal year is the time taken for the Sun to return to the same position with 
// respect to the stars of the celestial sphere. It is the orbital period of Earth, equal 
// to 365.25636042 mean solar days (31,558,149.540 seconds), that is 366.25636042 earth 
// rotations or sidereal days. The sidereal year is 20 minutes and 24 seconds longer than 
// the tropical year. The word "sidereal" means "relating to the stars". It derives from 
// the Latin sidus, meaning "star".

#define SECSYR 3.15581497632e7  


/*----- MASS OF THE SUN [Kg] -----*/

#define MASS_SUN 1.98844e30 


/*----- SPEED OF LIGHT [m/s] -----*/

#define LIGHT_SPEED 299792458.0


/*----- PARSEC [m] -----*/

#define PARSEC 3.08568025e16


/*----- Astronomical Unit [m] -----*/

#define AU 149597900000.0


/*----- Units and Physical Constants -----*/
  double CSPEED3;                                 // c^3
  
  double BHMASS_to_s;                     // Factor to convert ROK Times to Times in Seconds
  double BHMASS_to_m;                   // Factor to convert ROK Distances to Distances in Meters
  	
  //double BHMASS_to_s_o;                   // Auxiliar variable
  //double BHMASS_to_m_o;                 // Auxiliar variable


 // double GC_Conversion_from_meters_to_parsecs;          // Factor to convert lengths in meters to parsecs
  //double GC_Conversion_from_seconds_to_years;           // Factor to convert times in secons to sideral years
  

