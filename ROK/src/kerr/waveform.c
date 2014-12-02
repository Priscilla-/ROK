/*===============================================================================

    Computation of the Waveforms H_plus and H_cross up to Quadrupole order
 -------------------------------------------------------------------------------
                                              Created by Priscilla on 19/12/2012
 -------------------------------------------------------------------------------
 Last Update : 25.09.2013                                                            
 ===============================================================================*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/types.h>
#include<time.h>
#include <math.h>

#include "global_quantities_kerr.h"
#include "global_quantities_main.h"
#include "global_quantities_osc.h"
#include "physical_quantities.h"
#include "macros_2.h"


int analytic_r_costh(double*,double *,double *);
int checks_errors(int);
int dchi_dpsi_dphi_dt(double *, double, double, double *, double *,double *, double *);
int flat_space_XVA(double *x, double r, double costh, double dphidt, double drdt, double dcosthdt, double *Rxyz, double *Vxyz, double *Axyz);
int scalar_product(double *, double *, double *);
int vector_product(double *, double *, double *);
int vector_norm(double *, double *);


double *Realvector(const long, const long);
double **Realmatrix(const long, const long, const long, const long);

void free_Realvector(double *, const long, const long);
void free_Realmatrix(double **,const long , const long,const long , const long );

int waveform(double *H_p, double *H_x, double *x,double r, double costh)
{    
 
  FILE *data;
    
  double *Rxyz, *Vxyz, *Axyz;                    // Spatial location, velocity, and Axyzeleration (for waveform calculations)
  double vec2,Rxyz_dot_Axyz;                     // Some auxiliary quantities
  double **MQ;                                   // Components of the Mass Quadrupole
  double **Kronecker;                            // Kronecker's delta
  double dpsidt, dchidt,dphidt, dtaudt;
  double drdt, dcosthdt;
  long unsigned j,k;                             // Counters for Cartesian coordinates

/*----- Allocating memory for the Position, Velocity, and Acceleration of the Particle in BL associated Cartesian Coordinates -----*/
  Rxyz = Realvector(1,4); //Cartesian coordinates (X_p,Y_p,Z_p)
  Vxyz = Realvector(1,4); // velocity wrt the Boyer-Lindquist Time T_p
  Axyz = Realvector(1,4); // Acceleration wrt the BL time T_p

/*-- Computing analytical expressions for r(w_r) and cos(th)[w_th] and the radial and polar phases psi & chi --*/
 analytic_r_costh(x, &r, &costh);
    
/*----- Computing velocities d(Chi,Psi,Phi, dTau)/dt-----*/
 dchi_dpsi_dphi_dt(x, r, costh, &dpsidt, &dchidt,&dphidt, &dtaudt);
  
 drdt = dRdt;
 dcosthdt= dCOSTHdt;
    
 flat_space_XVA(x,r,costh,dphidt, drdt, dcosthdt, Rxyz, Vxyz, Axyz);

/*----- Scalar and vectorial products -----*/
  scalar_product(&vec2,Vxyz,Vxyz);
  scalar_product(&Rxyz_dot_Axyz,Rxyz,Axyz);
  
/*----- Computing the Mass Quadrupole -----*/
//----- Initializing the Kronecker's delta -----//
  Kronecker = Realmatrix(1, 3, 1, 3);
  for(j=1;j<=3;j++)
  {
    for(k=1;k<=3;k++)
    {
      if (k == j)
        Kronecker[j][k] = 1.0;
      else
        Kronecker[j][k] = 0.0;
    }
  }

  MQ = Realmatrix(1, 3, 1, 3);
  for(j=1;j<=3;j++)
  {
    for(k=1;k<=3;k++)
    {
	    if (j > k)
       MQ[j][k] = MQ[k][j]; // off-diagonal components
	    else
	     MQ[j][k] = MASS_SCO*( Rxyz[j]*Axyz[k] + Rxyz[k]*Axyz[j] + 2.0*Vxyz[j]*Vxyz[k]- 2.0*Kronecker[j][k]*(Rxyz_dot_Axyz + vec2)/3.0 );
    }
  }
    
/*-----------------------------------------------*/
/*----- Computing the Waveform Contribution -----*/
/*-----------------------------------------------*/

/*----- Initializing the Waveform Polarizations associated with this Observer -----*/
  *H_p = 0.0;
  *H_x = 0.0;
/*----- Adding the Mass Quadrupole -----*/
  for (j=1;j<=3;j++)
  {
    for (k=1;k<=3;k++)
    {
      *H_p += E_p[j][k]*MQ[j][k];
      *H_x += E_x[j][k]*MQ[j][k];
      
     // printf("%4,6e %4,6e %4,6e\n",E_p[j][k],E_x[j][k],MQ[j][k] );
    }
  }
/*----- Rescaling with the Distant from the Source to the Observer -----*/
 
	*H_p /= DL;
	*H_x /= DL;
  
// hLISA(n_sampling, t, h_I, h_II, H_p, H_x);
  
 free_Realmatrix(Kronecker, 1, 3, 1, 3);
 free_Realmatrix(MQ, 1, 3, 1, 3);

 free_Realvector(Rxyz, 1, 4);
 free_Realvector(Vxyz, 1, 4);
 free_Realvector(Axyz, 1, 4);
  
/*----- This is the end -----*/
  return 0;  
}


