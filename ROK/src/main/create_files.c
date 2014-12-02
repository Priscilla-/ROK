/*=====================================================================================
                                                                                       
 This routine creates the Directory Names, File Names, and Files that the Code uses.   
                                                                                       
--------------------------------------------------------------------------------------
                                                    Created by Priscilla on 19/12/2012
--------------------------------------------------------------------------------------
 Last Update: 19.01.13 by PCM
=====================================================================================*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"


int create_files(char *RUN_number)
{
  char shell_command[200];

/*----- Name Basis -----*/
  char run_dir[] = "data/Run_";
  
/*----- Name of the Run Directory -----*/
  strcpy(DirRun, run_dir);
  strcat(DirRun, RUN_number);
	
/*----- Creating Run Directory -----*/
  strcpy(shell_command, "mkdir -p -m700 ");
  strcat(shell_command, DirRun);

  system(shell_command);

/*----- Creating file to store run parameters -----*/
  strcpy(FN_parameters, DirRun);
  strcat(FN_parameters, "/run_parameters.txt");

/*----- Creating File to store  evoution related quantities -----*/
  strcpy(FN_evolution_AA, DirRun);
  strcat(FN_evolution_AA, "/Evolution_AA.txt");

  strcpy(FN_evolution_RR, DirRun);
  strcat(FN_evolution_RR, "/Evolution_RR.txt");
  
  
  strcpy(FN_fluxes, DirRun);
  strcat(FN_fluxes, "/Fluxes.txt");

/*---Creating file to store the orbital parameters and frequencies---*/
  strcpy(FN_frequencies, DirRun);
  strcat(FN_frequencies, "/Frequencies.txt");
   
/*---Creating file to store the orbital frequencies on Resonance---*/
  strcpy(FN_Res_frequencies, DirRun);
  strcat(FN_Res_frequencies, "/Frequencies_on_Resonance.txt");
    
/*---Creating file to store the GW waveforms---*/
  strcpy(FN_waveforms, DirRun);
  strcat(FN_waveforms, "/Waveforms.txt");

/*----- This is the end -----*/
  return 0;

}



