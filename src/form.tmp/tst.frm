*===This file was automatically generated in Python==*

off statistics;

*====includes===


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
CF F,J;

* == Gamma matrices == 
F G;
CF g;

* == Symbols ==
S g0, rh, rw, m0, a;  * // physical constants
S pow, pow2, sf, DIM; * // powers/scale factors
S db, dq, aa, bb;     * // bosonic/fermionic denominators
S n,nx,ny,nz,nt;      * // powers of k - used for integrating
S OP;                 * // momentum space operator
S iter,MAX,i;

CF cf;
autodeclare S c;
autodeclare S pw;


* == The expression==
L expr = 0+0+1/8*F(1,2,0,0,0,3)*rw^2*OP+3/8*F(1,2,0,0,1,2)*rw^2*OP+1/16*F(1,1,0,0,0,2)*rw^2*OP+3/16*F(1,1,0,0,1,1)*rw^2*OP+1/16*F(1,1,0,0,0,2)*rw^2*OP+3/16*F(1,1,0,0,1,1)*rw^2*OP+0+1/2*F(1,2,0,0,1,1)*OP+1/2*F(1,2,0,0,0,2)*OP+F(1,2,0,0,0,2)*OP+1/2*F(1,1,0,0,0,1)*OP-(0)-(0)-(1/2*F(1,2,0,0,0,2)*rw^2*OP+3/2*F(1,2,0,0,1,1)*rw^2*OP)-(F(1,1,0,0,0,1)*rw^2*OP)-(1/8*F(1,2,0,0,1,2)*OP)-(0)-(2*B(0,3,0,0,0,1)*OP+2*J(1,2,0,0,0,1)*OP)-(1/8*F(1,2,0,0,0,3)*OP)-(1/8*F(1,1,0,0,0,2)*OP)-(2*B(0,3,0,0,0,1)*OP+2*J(1,2,0,0,0,1)*OP);
id a^pow?!{-20,0}=0;

id F(pw1?,pw2?,pw3?,pw4?,pw5?,pw6?) = cf(pw1+pw2)*F(pw1,pw2,pw3,pw4,pw5,pw6);

contract;
print expr;
***********
*ANSWER
.end