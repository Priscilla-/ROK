/*==============================================================================
 
  Here we attempt to calculate the first few Fourier components of the forcing
  term for the evolution of the integrals of the motion.

  This is currently work in progress.
 
 -------------------------------------------------------------------------------
                                          Created by Christopher on 16/09/2013
 -------------------------------------------------------------------------------
 Last Update: 19/09/2013
 =============================================================================*/

#include<stdio.h>
#include<stdlib.h>

#include "global_quantities_main.h"
#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"
#include "macros_2.h"
#include "physical_quantities.h"

/*----Called subfunctions----*/
int PN_ELzQdot(double *, double *, double *, double *, double *, double , double , double , double , double , double );
int dchi_dpsi_dphi_dt(double *, double *, double , double , double *, double *, double *, double *);

/*===========================================================================
  Inputs are:
   x_Res - the vector of the resonant orbital parameters
   dxdt_dummy - this is passed to subfunctions, but not actually used by them
   n_r - the radial resonance integer
   n_th - the polar resonance integer
   n_Fourier - the number of Fourier components to be calculated
  Outputs are:
   aE_cos - array of cosine harmonic coefficients for energy forcing function
	bE_sin - array of sine harmonic coefficients for energy forcing function
   aL_cos - as aE_cos but for angular momentum
   bL_sin - as bE_sin but for angular momentum
   aQ_cos - as aE_cos but for Carter constant
   bQ_sin - as bE_sin but for Carter constant
  ===========================================================================*/
void Fourier_SF(double *x_Res, double *dxdt_dummy, int n_r, int n_th, int n_Fourier, double *aE_cos, double *bE_sin, double *aL_cos, double *bL_sin, double *aQ_cos, double *bQ_sin)
{

   int n_loop = 100; // Number of sampling points for Fourier series
   int m_loop = 100; // Number of sampling points for averaging
   int i_loop, j_loop, i_Fourier, i_int; // Loop indices
   double fE_tmp, fL_tmp, fQ_tmp, fE[n_loop + 1], fL[n_loop + 1], fQ[n_loop + 1]; // Forcing terms
   double q_r, q_th; // angle variables
   double dEdt_loop, dLdt_loop, dQdt_loop, dxdt_loop, x_loop[7], r_loop, costh_loop, dpsidt_loop, dchidt_loop, dphidt_loop, dtaudt_loop; // Arguments for PN_ELzQdot
   double h_step; // Step size for trapezium rule integration
   double recip_PI = 1.0 / PI;
   double nq_i, nq_j; // Phases for Fourier components at points i and j = i + 1 in integration
   double cos_i, sin_i, cos_j, sin_j; // Cosine and sine of nq_i and nq_j
   double da_E, da_L, da_Q, db_E, db_L, db_Q; // Integration elements for Fourier components
   double aE_int, aL_int, aQ_int, bE_int, bL_int, bQ_int; // Integrals for Fourier elements

   /*----To calculate the Fourier components we must evaluate the forcing term for all values of q_res----*/
   /*----For each value of q_res, q_r (or q_th) could have any value, it is necessary to average over these----*/

   /*---Orbital parameters----*/
   x_loop[1] = x_Res[1]; // phi
   x_loop[2] = x_Res[2]; // tau (not used)
   x_loop[3] = x_Res[3]; // semilatus rectum
   x_loop[4] = x_Res[4]; // eccentricity
   x_loop[5] = x_Res[5]; // inclination

   /*----Calculate forcing term over resonant orbit----*/
   for (i_loop = 0; i_loop < n_loop; i_loop++) 
   {
      /*----Average over auxiliary angle variable----*/
      fE_tmp = 0.0; // Running total for energy
      fL_tmp = 0.0; // Running total for angular momentum
      fQ_tmp = 0.0; // Running total for Carter constant

      for (j_loop = 0; j_loop < m_loop; j_loop++) // We pick q_r to average over
      {
         /*----Set up current position and velocity----*/
         /*----Calculate angle variables----*/
         q_r = 2.0 * PI * n_th * j_loop / m_loop; // q_r (we average over n_th cycles to ensure an even covering)
         q_th = (2.0 * PI * (i_loop / n_loop) - n_r * q_r) / n_th; // q_th (this should complete n_r cycles)

         x_loop[6] = q_r - floor(q_r / (2.0 * PI)); // radial angle variable in the range 0--2pi
         x_loop[7] = q_th - floor(q_th / (2.0 * PI)); // theta angle variable in the range 0--2pi

         analytic_r_costh(x_loop, &r_loop, &costh_loop); // Convert angle variables to BL coordinates

         dchi_dpsi_dphi_dt(dxdt_dummy, x_loop, r_loop, costh_loop, &dpsidt_loop, &dchidt_loop, &dphidt_loop, &dtaudt_loop); // Obtain velocities

         /*----Evaluate forcing terms at current position----*/
         /*----Use PN self-force expressions----*/
         PN_ELzQdot(&dEdt_loop, &dLdt_loop, &dQdt_loop, dxdt_dummy, x_loop, r_loop, costh_loop, dpsidt_loop, dchidt_loop, dphidt_loop, dtaudt_loop);

         fE_tmp += dEdt_loop;
         fL_tmp += dLdt_loop;
         fQ_tmp += dQdt_loop;

      }

      /*----Construct simple averages----*/   
      fE[i_loop] = fE_tmp / m_loop;
      fL[i_loop] = fL_tmp / m_loop;
      fQ[i_loop] = fQ_tmp / m_loop;
   }

	/*----Use periodicity----*/
   fE[n_loop] = fE[0];
   fL[n_loop] = fL[0];
   fQ[n_loop] = fQ[0];

   /*----Create Fourier series----*/
   /*----We do not calculate the 0th term as this gives adiabatic evolution----*/
   h_step = 2.0 * PI / n_loop;
   for (i_Fourier = 0; i_Fourier < n_Fourier; i_Fourier++)
   {
      /*----Integrate using trapezium rule----*/
      aE_int = da_E;
      bE_int = db_E;

      aL_int = da_L;
      bL_int = db_L;

      aQ_int = da_Q;
      bQ_int = db_Q;

      /*----Step through integration----*/
      for (i_int = 0; i_int < n_loop; i_int++)
      {
         /*----Terms for Fourier series----*/
         nq_i = (i_Fourier + 1) * 2.0 * PI * i_int;
         nq_j = (i_Fourier + 1) * 2.0 * PI * (i_int + 1);

         cos_i = cos(nq_i);
         cos_j = cos(nq_j);
         sin_i = sin(nq_i);
         sin_j = sin(nq_j);

         /*----Differential increments----*/
         da_E = 0.5 * h_step * (cos_i * fE[i_int] + cos_j * fE[i_int + 1]);
         db_E = 0.5 * h_step * (sin_i * fE[i_int] + sin_j * fE[i_int + 1]);

         da_L = 0.5 * h_step * (cos_i * fL[i_int] + cos_j * fL[i_int + 1]);
         db_L = 0.5 * h_step * (sin_i * fL[i_int] + sin_j * fL[i_int + 1]);

         da_Q = 0.5 * h_step * (cos_i * fQ[i_int] + cos_j * fQ[i_int + 1]);
         db_Q = 0.5 * h_step * (sin_i * fQ[i_int] + sin_j * fQ[i_int + 1]);

         /*----Integrals----*/
         aE_int += da_E;
         bE_int += db_E;

         aL_int += da_L;
         bL_int += db_L;

         aQ_int += da_Q;
         bQ_int += db_Q;

      }

      /*----Fourier components----*/
      aE_cos[i_Fourier] = recip_PI * aE_int;
      bE_sin[i_Fourier] = recip_PI * bE_int;

      aL_cos[i_Fourier] = recip_PI * aL_int;
      bL_sin[i_Fourier] = recip_PI * bL_int;

      aQ_cos[i_Fourier] = recip_PI * aQ_int;
      bQ_sin[i_Fourier] = recip_PI * bQ_int;
   }

}
