/*=====================================================================================
                                                                                    
  This Routine checks consistency and compatibility of different parts of the Code.    
                                                                                       
---------------------------------------------------------------------------------------
  "type" already used: 1,6
  "type" left free: 2-5,7-10
---------------------------------------------------------------------------------------
                                                    Created by Priscilla on 28/12/2012
---------------------------------------------------------------------------------------
  Last Update : 28.12.11
======================================================================================*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include "global_quantities_main.h"

int checks_main_errors(int type)
{

/*----- Switch that selects the particular check of a given call -----*/
  switch(type)
  {
    case 1 :
             /*----- Shows how to Run the Compiled ROK Code -----*/
	         printf("  |--------------------------------------------------------------------------------\n  | \n");
             printf("  | Usage: ./exe/rok Run_Number\n  | \n");
             printf("  | Format of Run_Number: XXXX (Four Digits!)\n  | \n");
             printf("  |--------------------------------------------------------------------------------\n\n");
             exit(1);
	        
    case 2 :
             /*-----  -----*/
	         printf("  |--------------------------------------------------------------------------------\n  | \n");
	         printf("  |--------------------------------------------------------------------------------\n\n");
             exit(1);
	        
    case 3 :
             /*-----  -----*/
	         printf("  |--------------------------------------------------------------------------------\n  |\n");
	         printf("  |--------------------------------------------------------------------------------\n");
             exit(1);
	        
    case 4 :
             /*-----  -----*/
	         printf("  |--------------------------------------------------------------------------------\n  | \n");
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
	         printf("  |--------------------------------------------------------------------------------\n  | \n");
             printf("  |--------------------------------------------------------------------------------\n\n");
             exit(1);
			 
	case 9:
	         /*-----  -----*/
			 printf("  |--------------------------------------------------------------------------------\n  |\n");
             printf("  |--------------------------------------------------------------------------------\n\n");
             exit(1);
			 	        
	case 10:
	         /*-----  -----*/
			 printf("  |--------------------------------------------------------------------------------\n  |\n");
             printf("  |--------------------------------------------------------------------------------\n\n");
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

