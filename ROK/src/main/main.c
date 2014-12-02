/*==============================================================================
    This is the ROK (Resonant Orbits in Kerr) code, its main goal is to search
    and analyze resonant geodesics in Kerr spacetime.
 -------------------------------------------------------------------------------
                                                            Priscilla Canizares
 -------------------------------------------------------------------------------
 Last Update : 06.10.13
 =============================================================================*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"
#include "macros_2.h"


int checks_main_errors(int);
int conversion_factors(void);
int create_files(char *);
int evolution_AA(FILE *,FILE *,FILE *,FILE *, double *x);
int evolution_RR(FILE *);
int find_resonances(void);
int initialization(void);
int modify_computational_parameters(void);
int reading_parameters(void);
int save_run_parameters(void);
int wrap_up();

double *Realvector(const long nl, const long nh);
void free_Realvector(double *v, const long nl, const long nh);

int main(int argc, char *argv[])
{
  double *x;
  char *run_number;
  FILE *data, *dataf, *dataflux, *dataGW;
  
/*----- Checking the number of Variables passed to the Executable [koala] -----*/
  if (argc != 2) checks_main_errors(1);
  if (strlen(argv[1]) != 4) checks_main_errors(1);
    
  run_number = argv[1];

/*----- Creating Run Directory and Data Files -----*/
  create_files(argv[1]);
    
/*----- Reading and storing parameters for this Run -----*/
  reading_parameters();
  
  save_run_parameters();

/*----- Compute Unit Conversion Factors -----*/
  conversion_factors();

/*----- Adjusting Time Variables and Counters -----*/
  modify_computational_parameters();
   
/*----- Searching Resonances -----*/
  if(FLAG_resonances == 'y')
   find_resonances();
  
/*----- Evolution using osculating formulation -----*/
  if(FLAG_OS=='y')
  {
    //---- Initialization of Global variables
    initialization();
    
    //--- variables to be evolved
   if(MAG_SF == 0)
    NVAR = 4;
   else        
    NVAR = 8;
   
   x = Realvector(1, NVAR);
      
   if(STORE_rate!= 0)
   {
      data   = fopen(FN_evolution_AA,"w");
      dataf  = fopen(FN_frequencies,"w");
      dataflux=fopen(FN_fluxes, "w");
      dataGW = fopen(FN_waveforms, "w");
   }
  
    evolution_AA(data, dataf, dataflux,dataGW, x);
    
   if(STORE_rate!= 0)
   {
      fclose(data);
      fclose(dataf);
      fclose(dataflux);
      fclose(dataGW);
   }
    
   free_Realvector(x, 1, NVAR);
  }
/*---- Evolution using Radiation Reaction effects----*/
  else
  {
//----Global variables initialization
    initialization();
    
    if(STORE_rate!= 0)
      data = fopen(FN_evolution_RR,"a");

    evolution_RR(data);

    if(STORE_rate!= 0)
    fclose(data);
  }
  
/*----- Wrap up -----*/
  wrap_up();
    
/*----- This is the end -----*/
  return 0;
}














