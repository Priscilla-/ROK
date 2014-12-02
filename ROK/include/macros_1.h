/* SOME MACROS */

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b), (dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
#define FMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b), (dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
#define LMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b), (dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
#define IMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b), (dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2)) 

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b), (dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))
#define FMIN(a,b) (dminarg1=(a),dminarg2=(b), (dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))
#define LMIN(a,b) (dminarg1=(a),dminarg2=(b), (dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))
#define IMIN(a,b) (dminarg1=(a),dminarg2=(b), (dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))


