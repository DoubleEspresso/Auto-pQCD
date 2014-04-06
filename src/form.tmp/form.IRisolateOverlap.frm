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
L expr = 256*sum(ix32)*g(ix32)*kh(4)^4*p(ix32)*q(ix32)^2*i_*db^-2*dq^-2*aa^2*bb^2*rh^2;


#call expandSumsOverlap()
.sort

#call irIsolateOverlap()
.sort 

#call projectOps()


contract;
print expr;
***********
*ANSWER
.end