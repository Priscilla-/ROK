/*****************************************************************************************/
/*                                                                                       */
/*          Computation and Output of the Waveforms H_plus and H_cross and of the        */
/*                                LISA response functions                                */
/* 					                                               						 */
/*****************************************************************************************/
/*                                                                                       */
/*  Last Update : 25.01.2013                                                             */
/*                                                                                       */
/*****************************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h> 
#include<sys/types.h>
#include<time.h>

#include "physical_quantities.h"
#include "global_quantities_koala.h"
#include "global_quantities_main.h"

int checks_errors(int);
// These functions are defined at the boton of the file
int scalar_product(double *, double *, double *);
int vector_product(double *, double *, double *);
int vector_norm(double *, double *);


int waveform(long unsigned n_sampling, double t_LISA_in_sec, double *h_I, double *h_II)
{    
 
  FILE *data;
    
  double pos[4], vec[4], acc[4];                 // Spatial location, velocity, and acceleration (for waveform calculations)
  double dacdt[4];                               // Time derivative of the acceleration)
  double pos_cross_vec[4], pos_cross_acc[4];     // Some auxiliary vectors
  double pos_cross_dacdt[4], vec_cross_acc[4];   //   "      "        "
  double gmf;                                    //Gamma Factor
  double pos2,vec2,pos_dot_vec,pos_dot_acc;      // Some auxiliary quantities
  double vec_dot_acc,pos_dot_dacdt;              //  "       "          "
    
  long unsigned j,k;                             // Counters for Cartesian coordinates
 
/**----- Cartesian Coordinate Systems -----**/
    struct vector3D
    {
        double x;
        double y;
        double z;
    };
    
    struct vector3D bl_r;              // Euclidean Cartesian coordinates (X_p,Y_p,Z_p).  These are Cartesian coordinates associated with Boyer-Lindquist coordinates in a trivial way 
    struct vector3D bl_v;              // Velocities of the Euclidean Cartesian coordinates (X_p,Y_p,Z_p) [see above] with respect to the Boyer-Lindquist Time T_p
    struct vector3D bl_ac;              // Accelerations of the Euclidean Cartesian coordinates (X_p,Y_p,Z_p) [see above] with respect to the Boyer-Lindquist Time T_p
    
    struct vector3D * bl_hd;           // Higher Derivatives of the Euclidean Cartesian coordinates (X_p,Y_p,Z_p) [see above] with respect to the Boyer-Lindquist Time T_p
/*----- Cartesian Position, Velocity, and Acceleration of the Particle in BL associated Cartesian Coordinates -----*/
  pos[1] = bl_r.x;   pos[2] = bl_r.y;  pos[3] = bl_r.z;
  vec[1] = bl_v.x;   vec[2] = bl_v.y;  vec[3] = bl_v.z;
  acc[1] = bl_ac.x;  acc[2] = bl_ac.y; acc[3] = bl_ac.z;

	
/*----- Higher-order time derivatives of the Cartesian Particle Position -----*/
  dacdt[1] = bl_hd[0].x;  dacdt[2] = bl_hd[0].y;  dacdt[3] = bl_hd[0].z;
	 	   
/*----- Some Derived Quantities -----*/
  scalar_product(&pos2,pos,pos);
  scalar_product(&vec2,vec,vec);
  scalar_product(&pos_dot_vec,pos,vec);	  
  scalar_product(&pos_dot_acc,pos,acc);
  scalar_product(&pos_dot_dacdt,pos,dacdt);
  scalar_product(&vec_dot_acc,vec,acc);
	
  vector_product(pos_cross_vec,pos,vec);
  vector_product(pos_cross_acc,pos,acc);
  vector_product(pos_cross_dacdt,pos,dacdt);
  vector_product(vec_cross_acc,vec,acc);


/*----- Gamma Factor -----*/
  gmf = 1.0/sqrt(1.0 - vec2);


/**----- Computing the Mass Quadrupole -----**/
  for(j=1;j<=3;j++)
  {
    for(k=1;k<=3;k++)
    {
	  if (j > k)
	    GV_Waveform.mass_quad[j][k] = GV_Waveform.mass_quad[k][j];
	  else
	    GV_Waveform.mass_quad[j][k] = GV_Mass_SCO*( pos[j]*acc[k] + pos[k]*acc[j] + 2.0*vec[j]*vec[k]- 2.0*GV_Kronecker[j][k]*(pos_dot_acc + vec2)/3.0 );
      	 
    }
  }
									

  //if ( N_higher_derivatives >= 1 )
  //{
    //checks_errors(10);
/**----- Computing the Mass Octupole -----**/
//      for (j=1;j<=3;j++)
//        for (k=1;k<=3;k++)
//          for (l=1;l<=3;l++)
//            GV_Waveform.mass_octo[j][k][l] = GV_Mass_SCO*( pos[j]*( dacdt[k]*pos[l] + pos[k]*dacdt[l] + 3.0*(acc[k]*vec[l]
//		                                   + vec[k]*acc[l]) ) + 3.0*vec[j]*( acc[k]*pos[l] + pos[k]*acc[l] + 2.0*vec[k]*vec[l] )
//		                                   + 3.0*acc[j]*( vec[k]*pos[l] + pos[k]*vec[l] ) + dacdt[j]*pos[k]*pos[l] 
//                                           - 0.2*( GV_Kronecker[j][k]*dacdt[l] + GV_Kronecker[j][l]*dacdt[k] + GV_Kronecker[k][l]*dacdt[j] )*pos2
//		                                   - 1.2*( GV_Kronecker[j][k]*acc[l] + GV_Kronecker[j][l]*acc[k] + GV_Kronecker[k][l]*acc[j] )*pos_dot_vec
//                                           - 1.2*( GV_Kronecker[j][k]*vec[l] + GV_Kronecker[j][l]*vec[k] + GV_Kronecker[k][l]*vec[j] )*( vec2 + pos_dot_acc )
//                                           - 0.4*( GV_Kronecker[j][k]*pos[l] + GV_Kronecker[j][l]*pos[k] + GV_Kronecker[k][l]*pos[j] )*( 3.0*vec_dot_acc + pos_dot_dacdt ) );
										  
/**----- Computing the Current Quadrupole -----**/
//      for (j=1;j<=3;j++)
//        for (k=1;k<=3;k++)
//          GV_Waveform.current_quad[j][k] = 0.5*GV_Mass_SCO*( acc[j]*pos_cross_vec[k] + pos_cross_vec[j]*acc[k] 
//	                                     + 2.0*( vec[j]*pos_cross_acc[k] + pos_cross_acc[j]*vec[k] )
//			                             + pos[j]*( vec_cross_acc[k] + pos_cross_dacdt[k] ) + ( vec_cross_acc[j] + pos_cross_dacdt[j] )*pos[k] );
 //}
	  


/*-----------------------------------------------*/
/*----- Computing the Waveform Contribution -----*/
/*-----------------------------------------------*/

/*----- Initializing the Waveform Polarizations associated with this Observer -----*/
  GV_LISA.h_plus = 0.0; 
  GV_LISA.h_cross = 0.0;  

  
/*----- Adding the Mass Quadrupole -----*/
  for (j=1;j<=3;j++)
  {
    for (k=1;k<=3;k++)
    {
      GV_LISA.h_plus += GV_LISA.polarization_tensor_plus[j][k]*GV_Waveform.mass_quad[j][k];
      GV_LISA.h_cross += GV_LISA.polarization_tensor_cross[j][k]*GV_Waveform.mass_quad[j][k];
    }
  }
    
    
/*----- Adding Multipoles with Higher-order Time Derivatives -----*/
//  if ( GV_N_higher_derivatives >= 1 )
  //{
    //checks_errors(10);
/**----- Adding the Current Quadrupole -----**/
//      for (i=1;i<=3;i++)
//      {
//        for (j=1;j<=3;j++)
//        {
//          for (k=1;k<=3;k++)
//          {
//            for (l=1;l<=3;l++)
//            {
//              GV_LISA.h_plus += 4.0*GV_LISA.polarization_tensor_plus[i][j]*GV_Levi_Civita[i][k][l]*GV_Waveform.current_quad[j][k]*GV_LISA.n_src[l]/3.0;
//              GV_LISA.h_cross += 4.0*GV_LISA.polarization_tensor_cross[i][j]*GV_Levi_Civita[i][k][l]*GV_Waveform.current_quad[j][k]*GV_LISA.n_src[l]/3.0;
//            }
//          }
//        }
//      }
	  
/**----- Adding the Mass Octupole -----**/
//      for (i=1;i<=3;i++)
//      {
//        for (j=1;j<=3;j++)
//        {
//           for (k=1;k<=3;k++)
//           {
//             GV_LISA.h_plus += GV_LISA.polarization_tensor_plus[i][j]*GV_Waveform.mass_octo[i][j][k]*GV_LISA.n_src[k]/3.0;
//             GV_LISA.h_cross += GV_LISA.polarization_tensor_cross[i][j]*GV_Waveform.mass_octo[i][j][k]*GV_LISA.n_src[k]/3.0;
//           }
//        }
//      }
 // }
			    
/**----- Rescaling with the Distant from the Source to the Observer -----**/
 
	GV_LISA.h_plus /= GV_LISA.d_l;
	GV_LISA.h_cross /= GV_LISA.d_l;
  
  
              
/*---- Creates a file where to store the GW polarizations h+ and hx ---*/                                                    
/*    if ( GV_Flag_waveforms == 'y')
    {        
		strcpy(fn_wpol, GV_Dn_run);
		strcat(fn_wpol, "/t_hpx_phase.txt");
        
        data_polarizations = fopen(fn_wpol,"a");
         fprintf(data_polarizations," %4.6e  %4.6e\n",t_LISA_in_sec, 0.5*PI + atan(fabs(GV_LISA.h_plus/GV_LISA.h_cross)) );
        fclose(data_polarizations);        
    }
*/    
    
    
/**----- Computing LISA/eLISA Response Functions and Saving WAVEFORMS and Response Functions -----**/
  if ( (GV_Integer_Computational_Parameter[11] == 1) || (GV_Integer_Computational_Parameter[11] == 4) )
  {   
    hLISA(n_sampling, t_LISA_in_sec, h_I, h_II, GV_LISA.h_plus, GV_LISA.h_cross);
  }
/**----- Computing aLIGO/ET Response Functions and Saving WAVEFORMS and Response Functions -----**/
  else if ( (GV_Integer_Computational_Parameter[11] == 2) || (GV_Integer_Computational_Parameter[11] == 3) )
  {
    h_ground_based(n_sampling, t_LISA_in_sec, h_I, h_II, GV_LISA.h_plus, GV_LISA.h_cross);
  }


/*----- This Routine ends Here! -----*/
  return 0;  
}



/*===== SOME AUXILIARY ROUTINES =====*/

/**===== SCALAR PRODUCT =====**/
int scalar_product(double *result, double *vector_1, double *vector_2)
{
  *result = vector_1[1]*vector_2[1] + vector_1[2]*vector_2[2] + vector_1[3]*vector_2[3];
  
  return 0;
}

/**===== VECTOR PRODUCT =====**/
int vector_product(double *result, double *vector_1, double *vector_2)
{
  int j,k;

  result[1] = 0.0;  result[2] = 0.0;  result[3] = 0.0;
  
  for (j=1;j<=3;j++)
  {
    for (k=1;k<=3;k++)
	{
	  result[1] += GV_Levi_Civita[1][j][k]*(vector_1[j])*(vector_2[k]);
	  result[2] += GV_Levi_Civita[2][j][k]*(vector_1[j])*(vector_2[k]);
	  result[3] += GV_Levi_Civita[3][j][k]*(vector_1[j])*(vector_2[k]);
	}
  }
 
  return 0;
}

/**===== VECTOR NORM =====**/
int vector_norm(double *result, double *vector)
{
  
  *result = sqrt( vector[1]*vector[1] + vector[2]*vector[2] + vector[3]*vector[3] );
  
  return 0;
}