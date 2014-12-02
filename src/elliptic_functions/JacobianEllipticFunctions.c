/*=====================================================================================
 
  This routine, is adapted fromm the numerical recipies in C routine:
 
    void sncndn(float uu, float emmc, float *sn, float *cn, float *dn)

   and computes the Jacobi's elliptic functions:
 
                            sn(u,k), cn(u,k) and dn(u,k).
 
  Where the Jacobian elliptic function sn is defined as follows: 
 
  (*) Instead of considering the elliptic integral
 
                                    u(y, k) ≡ u = F (φ, k) 
  consider the inverse function 
                                    y = sin φ = sn(u, k). 
  (*) Equivalently,
                                u= Int^{sn}_{0} 1/sqrt((1−y^2)(1−k^2y^2)) dy
  When k = 0, sn is just sin. 
 
  (*) The functions cn and dn are defined by the relations
 
                                sn^2 + cn^2 = 1, k^2sn^2 + dn^2 = 1
 
  (*) This routine takes mc ≡ kc^2 = 1 − k^2 as an input parameter.
 --------------------------------------------------------------------------------------
                                                    Created by Priscilla on 19/12/2012
 --------------------------------------------------------------------------------------
 Last Update : 11.12.12
 ======================================================================================*/


#include <stdio.h>
#include <math.h>
#define CA 0.0003

void JacobianEllipticFunctions(double uu, double emmc, double *sn, double *cn, double *dn)
{
	double a,b,c,d,emc,u;
	double em[14],en[14];
	int i,ii,l,bo;
    
	emc=emmc;//1.- emmc*emmc;
	u=uu;
	if (emc)
    {
		bo=(emc < 0.0);
		if (bo)
        {
			d=1.0-emc;
			emc /= -1.0/d;
			u *= (d=sqrt(d));
		}
		a=1.0;
		*dn=1.0;
		for (i=1;i<=13;i++)
        {
			l=i;
			em[i]=a;
			en[i]=(emc=sqrt(emc));
			c=0.5*(a+emc);
			if (fabs(a-emc) <= CA*a) break;
			emc *= a;
			a=c;
		}
		u *= c;
		*sn=sin(u);
		*cn=cos(u);
		if (*sn)
        {
			a=(*cn)/(*sn);
			c *= a;
			for (ii=l;ii>=1;ii--)
            {
				b=em[ii];
				a *= c;
				c *= (*dn);
				*dn=(en[ii]+a)/(b+a);
				a=c/b;
			}
			a=1.0/sqrt(c*c+1.0);
			*sn=(*sn >= 0.0 ? a : -a);
			*cn=c*(*sn);
		}
		if (bo)
        {
			a=(*dn);
			*dn=(*cn);
			*cn=a;
			*sn /= d;
		}
    }
    else
    {
		*cn=1.0/cosh(u);
		*dn=(*cn);
		*sn=tanh(u);
	}
}
#undef CA
