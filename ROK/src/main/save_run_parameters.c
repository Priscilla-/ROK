/*==============================================================================
 
Here we write in a file the parameters employed in the Fisher Matrix computations
 
-------------------------------------------------------------------------------
  Last Update : 01.04.12
===============================================================================*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"
#include "macros_2.h"



int save_run_parameters(void)
{
    FILE *data;
    
 data = fopen(FN_parameters,"a");
    fprintf(data,"======================================================================================\n\n");
    fprintf(data,"          SYSTEM PARAMETERS     \n\n");
    fprintf(data," Mass_BH = %4.6e   \n", EMRI_Parameter[1]);
    fprintf(data," Spin = %4.6e [M.] \n", EMRI_Parameter[2]);
    fprintf(data," Mu=m/M = %4.6e    \n", EMRI_Parameter[3]);
    fprintf(data," e_o = %4.6e\n",   EMRI_Parameter[4]);
    fprintf(data," p_o = %4.6e [M.]  \n",EMRI_Parameter[5]);
    fprintf(data," Theta_inc = %4.6e \n", EMRI_Parameter[6]);
    fprintf(data," Theta_Source = %4.6e \n", EMRI_Parameter[7]);
    fprintf(data," Phi_Source = %4.6e	  \n", EMRI_Parameter[8]);
    fprintf(data," Theta_Kerr = %4.6e	  \n", EMRI_Parameter[9]);
    fprintf(data," Phi_Kerr = %4.6e      \n", EMRI_Parameter[10]);
    fprintf(data," D_L/mu = %4.6e   \n", EMRI_Parameter[11]);
    fprintf(data," Psi_o = %4.6e	      \n", EMRI_Parameter[12]);
    fprintf(data," Chi_o = %4.6e    \n", EMRI_Parameter[13]);
    fprintf(data," PHI_o = %4.6e   \n", EMRI_Parameter[14]);
    fprintf(data,"----------------------------------------------------------------------------------------- \n");
    fprintf(data, " Searching Resonances [y/n]: %d \n",FLAG_resonances);
    fprintf(data, " k=%d,  n=%d \n", Resorder[0], Resorder[1]);
    fprintf(data, " e_f=%4.6e \n", EMRI_Parameter[15]);
    fprintf(data, " de =%4.6e\n", Real_Parameter[6]);
    fprintf(data, " p_f=%4.6e\n", EMRI_Parameter[16]);
    fprintf(data, " dp =%4.6e\n", Real_Parameter[7]);
    fprintf(data,"----------------------------------------------------------------------------------------- \n");
    fprintf(data," Include Radiation Reaction: %c   \n\n", FLAG_RR);
    fprintf(data," Evolution Time: T =  %1.2e   [y] \n\n", TIME_yrs);
    fprintf(data," Sampling time: t_sam = %4.4e  [s]\n\n", DTsec  );
    fprintf(data,"======================================================================================\n\n");
 fclose(data);
                    
/*----This is the end -----*/
    return 0;
}	
                    