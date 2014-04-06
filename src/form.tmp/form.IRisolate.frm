*===This file was automatically generated in Python==*

off statistics;

*====includes===
#include ../../inc/diracSimplify.h
#include ../../inc/symmetryTables.h
#include ../../inc/replaceBterms.h
#include ../../inc/projectOps.h
#include ../../inc/irIsolate.h
#include ../../inc/expandSums.h
#include ../../inc/expandSumsOverlap.h
#include ../../inc/formIntegrate.h
#include ../../inc/diracTrace.h

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
CF c;

* == Symbols ==
S g0, rh, rw, m0, a;  * // physical constants
S pow, pow2, sf, DIM; * // powers/scale factors
S db, dq, aa, bb;     * // bosonic/fermionic denominators
S n,nx,ny,nz,nt;      * // powers of k - used for integrating
S OP;                 * // momentum space operator

* == The expression==
L expr = 64*1*sum(ix26)*sum(ix28)*g(ix28)*sum_(ix34,1,4,kh( ix34)^2)*kh(4)^2*p(ix28)*q(ix26)^2*i_*rw^2*db^-2*dq^-1*aa^2*bb^2*rh^2;


#call expandSums()
.sort

#call irIsolate()
.sort 

#call projectOps()


contract;
print expr;
***********
*ANSWER
.end