/*=====================================================================================
 
  This routine computes is based on "Numerical Recipies in C" and computes:
 
  @ K(u) first kind elliptic function
  @ E(u) second kind elliptic function
  @ P(u) Third kind  elliptic function
  @ rj, rf, rd & rc: Carlson's elliptic integrals
 --------------------------------------------------------------------------------------
                                                    Created by Priscilla on 19/12/2012
 --------------------------------------------------------------------------------------
  Last Update : 11.12.12
 ======================================================================================*/

#include <stdio.h>
#include <math.h>
#define NRANSI

#include "global_quantities_kerr.h"

#include "macros_1.h"
#include "macros_2.h"

double rf(double, double, double);
double rd(double, double, double);
double rj(double, double, double, double);
double rc(double, double);

/*---Legendre elliptic integral of the 1st kind----*/ 
double ellK(double phi, double ak)
{
	double rf(double x, double y, double z);
	double s;
    
	s=sin(phi);
	return s*rf(SQR(cos(phi)),(1.0-s*ak)*(1.0+s*ak),1.0);
}
#undef NRANSI

/*---Legendre elliptic integral of the 2nd kind----*/ 
double ellE(double phi, double ak)
{
	double rd(double x, double y, double z);
	double rf(double x, double y, double z);
	double cc,q,s;
    
	s=sin(phi);
	cc=SQR(cos(phi));
	q=(1.0-s*ak)*(1.0+s*ak);
	return s*(rf(cc,q,1.0)-(SQR(s*ak))*rd(cc,q,1.0)/3.0);
}
#undef NRANSI

/*---Legendre elliptic integral of the 3rd kind----*/ 
double ellP(double phi, double en, double ak)
{
	double rf(double x, double y, double z);
	double rj(double x, double y, double z, double p);
	double cc,enss,q,s;
    
	s=sin(phi);
	enss=en*s*s;
	cc=SQR(cos(phi));
	q=(1.0-s*ak)*(1.0+s*ak);
	return s*(rf(cc,q,1.0)-enss*rj(cc,q,1.0,1.0+enss)/3.0);
}
#undef NRANSI

/*=================
  Auxiliar functions
 ==================*/

/*---Carlson’s elliptic integral of the third kind---*/
#define ERRTOL 0.05
#define TINY 2.5e-13
#define BIG 9.0e11
#define C1 (3.0/14.0)
#define C2 (1.0/3.0)
#define C3 (3.0/22.0)
#define C4 (3.0/26.0)
#define C5 (0.75*C3)
#define C6 (1.5*C4)
#define C7 (0.5*C2)
#define C8 (C3+C3)

double rj(double x, double y, double z, double p)
{
	double rc(double x, double y);
	double rf(double x, double y, double z);
	double a,alamb,alpha,ans,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,
    ed,ee,fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt;
    
	if (FMIN(FMIN(x,y),z) < 0.0 || FMIN(FMIN(x+y,x+z),FMIN(y+z,fabs(p))) < TINY
		|| FMAX(FMAX(x,y),FMAX(z,fabs(p))) > BIG) EVAL_arg=1;
   //   {printf("invalid arguments in elliptic integral rj\n");}
	sum=0.0;
	fac=1.0;
	if (p > 0.0) {
		xt=x;
		yt=y;
		zt=z;
		pt=p;
	} else {
		xt=FMIN(FMIN(x,y),z);
		zt=FMAX(FMAX(x,y),z);
		yt=x+y+z-xt-zt;
		a=1.0/(yt-p);
		b=a*(zt-yt)*(yt-xt);
		pt=yt+b;
		rho=xt*zt/yt;
		tau=p*pt/yt;
		rcx=rc(rho,tau);
	}
	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		alpha=SQR(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz);
		beta=pt*SQR(pt+alamb);
		sum += fac*rc(alpha,beta);
		fac=0.25*fac;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		pt=0.25*(pt+alamb);
		ave=0.2*(xt+yt+zt+pt+pt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
		delp=(ave-pt)/ave;
	} while (FMAX(FMAX(fabs(delx),fabs(dely)),
                  FMAX(fabs(delz),fabs(delp))) > ERRTOL);
	ea=delx*(dely+delz)+dely*delz;
	eb=delx*dely*delz;
	ec=delp*delp;
	ed=ea-3.0*ec;
	ee=eb+2.0*delp*(ea-ec);
	ans=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))
                     +delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave));
	if (p <= 0.0) ans=a*(b*ans+3.0*(rcx-rf(xt,yt,zt)));
	return ans;
}
#undef ERRTOL
#undef TINY
#undef BIG
#undef C1
#undef C2
#undef C3
#undef C4
#undef C5
#undef C6
#undef C7
#undef C8
#undef NRANSI

/*---Computes Carlson’s elliptic integral of the first kind---*/
#define ERRTOL 0.08
#define TINY 1.5e-38
#define BIG 3.0e37
#define THIRD (1.0/3.0)
#define C1 (1.0/24.0)
#define C2 0.1
#define C3 (3.0/44.0)
#define C4 (1.0/14.0)

double rf(double x, double y, double z)
{
	double alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt;
    
	if (FMIN(FMIN(x,y),z) < 0.0 || FMIN(FMIN(x+y,x+z),y+z) < TINY ||
		FMAX(FMAX(x,y),z) > BIG) EVAL_arg =1;
    //{printf("invalid arguments in elliptic integral rf\n"); }
	xt=x;
	yt=y;
	zt=z;
	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		ave=THIRD*(xt+yt+zt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
	} while (FMAX(FMAX(fabs(delx),fabs(dely)),fabs(delz)) > ERRTOL);
	e2=delx*dely-delz*delz;
	e3=delx*dely*delz;
	return (1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave);
}
#undef ERRTOL
#undef TINY
#undef BIG
#undef THIRD
#undef C1
#undef C2
#undef C3
#undef C4
#undef NRANSI


/*---Computes Carlson’s elliptic integral of the second kind---*/
#define ERRTOL 0.05
#define TINY 1.0e-25
#define BIG 4.5e21
#define C1 (3.0/14.0)
#define C2 (1.0/6.0)
#define C3 (9.0/22.0)
#define C4 (3.0/26.0)
#define C5 (0.25*C3)
#define C6 (1.5*C4)

double rd(double x, double y, double z)
{
	double alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,
    sqrtz,sum,xt,yt,zt;
    
	if (FMIN(x,y) < 0.0 || FMIN(x+y,z) < TINY || FMAX(FMAX(x,y),z) > BIG)
        EVAL_arg=1;
   // {printf("invalid arguments in elliptic integral rd\n");}
	xt=x;
	yt=y;
	zt=z;
	sum=0.0;
	fac=1.0;
	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		sum += fac/(sqrtz*(zt+alamb));
		fac=0.25*fac;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		ave=0.2*(xt+yt+3.0*zt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
	} while (FMAX(FMAX(fabs(delx),fabs(dely)),fabs(delz)) > ERRTOL);
	ea=delx*dely;
	eb=delz*delz;
	ec=ea-eb;
	ed=ea-6.0*eb;
	ee=ed+ec+ec;
	return 3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)
                        +delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave));
}
#undef ERRTOL
#undef TINY
#undef BIG
#undef C1
#undef C2
#undef C3
#undef C4
#undef C5
#undef C6
#undef NRANSI

/*---Computes Carlson’s degenerate elliptic integral---*/
#define ERRTOL 0.04
#define TINY 1.69e-38
#define SQRTNY 1.3e-19
#define BIG 3.e37
#define TNBG (TINY*BIG)
#define COMP1 (2.236/SQRTNY)
#define COMP2 (TNBG*TNBG/25.0)
#define THIRD (1.0/3.0)
#define C1 0.3
#define C2 (1.0/7.0)
#define C3 0.375
#define C4 (9.0/22.0)

double rc(double x, double y)
{
	double alamb,ave,s,w,xt,yt;
	if (x < 0.0 || y == 0.0 || (x+fabs(y)) < TINY || (x+fabs(y)) > BIG ||
		(y<-COMP1 && x > 0.0 && x < COMP2)) EVAL_arg=1;
   // {printf("invalid arguments in elliptic integral rc\n");}
	if (y > 0.0) {
		xt=x;
		yt=y;
		w=1.0;
	} else {
		xt=x-y;
		yt = -y;
		w=sqrt(x)/sqrt(xt);
	}
	do {
		alamb=2.0*sqrt(xt)*sqrt(yt)+yt;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		ave=THIRD*(xt+yt+yt);
		s=(yt-ave)/ave;
	} while (fabs(s) > ERRTOL);
	return w*(1.0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave);
}
#undef ERRTOL
#undef TINY
#undef SQRTNY
#undef BIG
#undef TNBG
#undef COMP1
#undef COMP2
#undef THIRD
#undef C1
#undef C2
#undef C3
#undef C4
#undef NRANSI






