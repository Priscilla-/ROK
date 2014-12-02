//
//  algebra.c
//  X_ROK_V1.7
//
//  Created by Priscilla on 17/08/2013.
//
//


#include<stdlib.h>
#include<stdio.h>
#include<string.h>


#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"

#define K_MAXIMUM 8               /*  Maximum row number in the extrapolation              */
#define KSEC (K_MAXIMUM+2)        /*                                                       */
#define SAFETY_1 0.25             /*  Safety factors                                       */
#define SAFETY_2 0.7              /*                                                       */
#define TINY 1.0e-30              /*  Prevents division by zero                            */
#define RED_MIN 0.7               /*  Minimum factor for step size reduction               */
#define RED_MAX 1.0e-5            /*  Maximum factor for step size reduction               */
#define SCALE_MAX 0.1             /*  1/SCALE_MAX = Maximum factor for step increasing     */

#define NR_END 1
#define FREE_ARG char*


//===== MATRIX & VECTORS from Scott's Code=====//
double *Realvector(const long nl, const long nh)
{
  double *v;
  
  v = (double *)malloc((size_t)((nh - nl + 1 + NR_END)*sizeof(double)));
  if (!v) printf("allocation failure in Realvector() \n");
  return v - nl + NR_END;
}

double **Realmatrix(const long nrl, const long nrh, const long ncl, const long nch)
{
  long i, nrow = nrh - nrl + 1, ncol= nch - ncl + 1;
  double **m;
  
  m = (double **)malloc((size_t)((nrow + NR_END)*sizeof(double*)));
  if (!m) printf("allocation failure 1 in Realmatrix()\n");
  m += NR_END;
  m -= nrl;
  
  m[nrl] = (double *)malloc((size_t)((nrow*ncol + NR_END)*sizeof(double)));
  if (!m[nrl]) printf("allocation failure 2 in Realmatrix()\n");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;
  
  return m;
}

void free_Realvector(double *v, const long nl,const long nh)
{
  if (nh < nl) printf("nh < nl in free_Realvector! \n");
  free((FREE_ARG)(v + nl - NR_END));
}

void free_Realmatrix(double **m,const long nrl, const long nrh,const long ncl, const long nch)
{
  if (nrh < nrl || nch < ncl)
    printf("Arguments out of whack in free_Realmatrix!");
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}
