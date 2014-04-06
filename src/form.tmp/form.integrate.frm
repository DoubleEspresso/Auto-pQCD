*===This file was automatically generated in Python==*

off statistics;

*====includes===
#include ../../inc/diracSimplify.h
#include ../../inc/symmetryTables.h
#include ../../inc/replaceBterms.h
#include ../../inc/projectOps.h
#include ../../inc/irIsolate.h
#include ../../inc/expandSums.h
#include ../../inc/formIntegrate.h

*==lorentz indices===
autodeclare index ix =0;
autodeclare index ex =0;

*==polarization indices===
autodeclare index sx = 0;

*==vectors==
CF p,p1,p2,p3,p4;  * // Fermion momenta
CF k;              * // Loop momenta 
CF q,q1,q2,q3,q4;  * // gauge field momenta 
CF d;              * // light-like vector 

* == Special functions == 
CF sum,SUM;         * // index summation helpers 
CF s,S,m,M,b,B,c,C; * // trig functions of loop-momenta

* == Gamma matrices == 
F G;
CF g;

* == Symbols ==
S g0, rw, m0, a;      * // physical constants
S pow, pow2, sf, DIM; * // powers/scale factors
S db, dq, aa, bb;     * // bosonic/fermionic denominators
S n,nx,ny,nz,nt;      * // powers of k - used for integrating
S OP;                 * // momentum space operator

* == The expression==
L expr = 1/8*F(2,2,0,0,0,4)*i_*OP+3/2*F(2,2,0,0,1,3)*i_*OP+9/8*F(2,2,0,0,2,2)*i_*OP+9/2*F(2,2,0,1,1,2)*i_*OP+3/4*F(2,2,1,1,1,1)*i_*OP+1/2*F(2,2,0,0,0,3)*i_*OP+9/2*F(2,2,0,0,1,2)*i_*OP+3*F(2,2,0,1,1,1)*i_*OP+1/8*F(2,1,0,0,0,3)*i_*OP+9/8*F(2,1,0,0,1,2)*i_*OP+3/4*F(2,1,0,1,1,1)*i_*OP+1/8*F(2,2,0,0,2,3)*i_*OP+1/8*F(2,2,0,1,2,2)*i_*OP+2*F(2,2,0,0,1,2)*i_*OP+2*F(2,2,0,1,1,1)*i_*OP+1/4*F(2,1,0,0,0,2)*i_*OP+3/4*F(2,1,0,0,1,1)*i_*OP+1/8*F(2,2,0,0,2,3)*i_*OP+1/8*F(2,2,0,1,2,2)*i_*OP+2*F(2,2,0,0,1,2)*i_*OP+2*F(2,2,0,1,1,1)*i_*OP+1/4*F(2,1,0,0,0,2)*i_*OP+3/4*F(2,1,0,0,1,1)*i_*OP+1/32*F(2,2,0,0,0,7)*i_*OP+15/32*F(2,2,0,0,1,6)*i_*OP+33/32*F(2,2,0,0,2,5)*i_*OP+45/32*F(2,2,0,0,3,4)*i_*OP+15/8*F(2,2,0,1,1,5)*i_*OP+105/16*F(2,2,0,1,2,4)*i_*OP+15/4*F(2,2,0,1,3,3)*i_*OP+75/16*F(2,2,0,2,2,3)*i_*OP+15/8*F(2,2,1,1,1,4)*i_*OP+15/2*F(2,2,1,1,2,3)*i_*OP+45/16*F(2,2,1,2,2,2)*i_*OP+1/2*F(2,2,0,0,0,3)*i_*OP+9/2*F(2,2,0,0,1,2)*i_*OP+3*F(2,2,0,1,1,1)*i_*OP+1/8*F(2,2,0,0,0,4)*i_*OP+3/8*F(2,2,0,0,2,2)*i_*OP+2*B(0,4,0,0,0,2)*i_*OP+6*B(0,4,0,0,1,1)*i_*OP+2*J(2,2,0,0,0,2)*i_*OP+6*J(2,2,0,0,1,1)*i_*OP+1/16*F(2,2,0,0,0,5)*i_*OP+3/16*F(2,2,0,0,1,4)*i_*OP+1/16*F(2,2,0,0,0,5)*i_*OP+3/16*F(2,2,0,0,1,4)*i_*OP+F(2,2,0,0,0,3)*i_*OP+3*F(2,2,0,0,1,2)*i_*OP+F(2,2,0,0,0,3)*i_*OP+3*F(2,2,0,0,1,2)*i_*OP+1/4*F(2,2,0,0,0,4)*i_*OP+4*B(0,4,0,0,0,2)*i_*OP+4*J(2,2,0,0,0,2)*i_*OP+1/8*F(2,2,0,0,2,3)*i_*OP+1/8*F(2,2,0,1,2,2)*i_*OP+2*F(2,2,0,0,1,2)*i_*OP+2*F(2,2,0,1,1,1)*i_*OP+1/4*F(2,1,0,0,0,2)*i_*OP+3/4*F(2,1,0,0,1,1)*i_*OP+1/8*F(2,2,0,0,2,3)*i_*OP+1/8*F(2,2,0,1,2,2)*i_*OP+2*F(2,2,0,0,1,2)*i_*OP+2*F(2,2,0,1,1,1)*i_*OP+1/4*F(2,1,0,0,0,2)*i_*OP+3/4*F(2,1,0,0,1,1)*i_*OP+1/2*F(2,2,0,0,2,2)*i_*OP+1/16*F(2,1,0,0,2,2)*i_*OP+8*B(0,4,0,0,1,1)*i_*OP+8*J(2,2,0,0,1,1)*i_*OP+F(2,1,0,0,1,1)*i_*OP+1/16*F(2,1,0,0,0,4)*i_*OP+1/2*F(2,2,0,0,2,2)*i_*OP+1/16*F(2,1,0,0,2,2)*i_*OP+F(2,1,0,0,0,2)*i_*OP+8*B(0,4,0,0,1,1)*i_*OP+8*J(2,2,0,0,1,1)*i_*OP+F(2,1,0,0,1,1)*i_*OP;

 

#call formIntegrate(expr)
.sort


contract;
print expr;
***********
*ANSWER
.end