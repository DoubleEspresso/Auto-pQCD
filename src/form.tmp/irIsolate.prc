#procedure irIsolate()

*// ======================================================================
*// this routine iteratively separates all IR-divergent
*// integrals and shifts them into newly defined B integrals.
*// this is done by 'splitting' the quark propagators
*// 
*//  1/q*f(q) = (1/b + (1/q-1/b)) * f(q)
*// 
*// collect the 1/q-1/b piece into a new J-integral, and label
*// the 1/b*f(q) piece by B.  All B-factors are IR-divergent and
*// can be computed analytically.  The J-integrals are IR-finite 
*// after enough iterations, and can be computed numerically.
*//
*// an integral is IR-divergent if it fails to satisfy the counting:
*// 
*// sum of denominator (q+b) powers - sum of numerator kh-powers < 2
*//
*// one iteration will reduce the IR-divergence by one.  In general, it may
*// be necessary to iterate more than once for some integrals.
*//
*//
*//====================================================================


*// we start by scaling all integrals by an factor which will check for IR divergences
*// the number of iterations is defined to be, 
*// iter = d_sum - (kh_sum + 2) + 1, cases with more than one iteration will need
*// to be handled.

*// to simplify the FORM coding, we have only coded reduction for one iteration
*// Integrals needing more than on iteration will have IRCHECK factors in front 
*// of them anyway.  

id F(qPow?, bPow?, pw1?, pw2?, pw3?, pw4?) = IRCHECK(bPow+qPow-pw1-pw2-pw3-pw4-1)
                                             *F(qPow, bPow, pw1, pw2, pw3, pw4);  

*=========================================================================
*// now we make conditional replacements based on the power of IRCHECK

*// case for one iteration - the splitting here is a little different
*// for overlap fermions, but can be handled in mathematica integration
*// routines.

id IRCHECK(1) * F(qPow?,bPow?,pw1?,pw2?,pw3?,pw4?) = 
     B(0,bPow+qPow,pw1,pw2,pw3,pw4) + J(qPow,bPow,pw1,pw2,pw3,pw4);

* // now re-label all IRCHECKS that don't require reduction.
id IRCHECK(0)  = 1; id IRCHECK(-1) = 1; id IRCHECK(-2) = 1; id IRCHECK(-3) = 1;
id IRCHECK(-4) = 1; id IRCHECK(-5) = 1; id IRCHECK(-6) = 1; id IRCHECK(-7) = 1; 
id IRCHECK(-8) = 1; id IRCHECK(-9) = 1; id IRCHECK(-10) = 1;id IRCHECK(-11) = 1;
id IRCHECK(-12) = 1; id IRCHECK(-13) = 1; id IRCHECK(-14) = 1;id IRCHECK(-15) = 1;
id IRCHECK(-16) = 1; id IRCHECK(-17) = 1; id IRCHECK(-18) = 1;id IRCHECK(-19) = 1;
id IRCHECK(-20) = 1; id IRCHECK(-21) = 1;
id IRCHECK(-22) = 1; id IRCHECK(-23) = 1; id IRCHECK(-24) = 1;id IRCHECK(-25) = 1;
id IRCHECK(-26) = 1; id IRCHECK(-27) = 1; id IRCHECK(-28) = 1;id IRCHECK(-29) = 1;
id IRCHECK(-30) = 1; id IRCHECK(-31) = 1; id IRCHECK(-32) = 1;id IRCHECK(-33) = 1;

