//
//  algebra.c
//  X_ROK_V1.7
//
//  Created by Priscilla on 17/08/2013.
//
//  modified 15/09/13


#include<stdlib.h>
#include<stdio.h>
#include<string.h>


#include "global_quantities_kerr.h"
#include "global_quantities_osc.h"
#include "global_quantities_main.h"

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

//===== SCALAR PRODUCT =====//
int scalar_product(double *result, double *vector_1, double *vector_2)
{
  *result = vector_1[1]*vector_2[1] + vector_1[2]*vector_2[2] + vector_1[3]*vector_2[3];
  
  return 0;
}

//===== VECTOR PRODUCT =====//
int vector_product(double *result, double *vector_1, double *vector_2)
{
  int i,j,k;
  
  /*----- Initialization of the 3-D Levi-Civita antisymmetric Symbol -----*/
  for (i=1;i<=3;i++)
    for (j=1;j<=3;j++)
      for (k=1;k<=3;k++)
        Levi_Civita[i][j][k] = 0.0;
  
  Levi_Civita[1][2][3] = 1.0;  Levi_Civita[2][3][1] = 1.0;  Levi_Civita[3][1][2] = 1.0;
  Levi_Civita[3][2][1] = -1.0; Levi_Civita[2][1][3] = -1.0;  Levi_Civita[1][3][2] = -1.0;
  

  
  result[1] = 0.0;  result[2] = 0.0;  result[3] = 0.0;
  
  for (j=1;j<=3;j++)
  {
    for (k=1;k<=3;k++)
    {
      result[1] += Levi_Civita[1][j][k]*(vector_1[j])*(vector_2[k]);
      result[2] += Levi_Civita[2][j][k]*(vector_1[j])*(vector_2[k]);
      result[3] += Levi_Civita[3][j][k]*(vector_1[j])*(vector_2[k]);
    }
  }
  
  return 0;
}

//===== VECTOR NORM =====//
int vector_norm(double *result, double *vector)
{
  
  *result = sqrt( vector[1]*vector[1] + vector[2]*vector[2] + vector[3]*vector[3] );
  
  return 0;
}



