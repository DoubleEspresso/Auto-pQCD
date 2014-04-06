off statistics;
format mathematica;

*-------------------------------------------------
* the color indices
Autodeclare Index a,b,c,d,e,f,g,h;

*-------------------------------------------------
* the indices on T-functions
Autodeclare Index i,j,l,m,n,o,r,s,t,dum;

*-------------------------------------------------
* the lorentz indices
Autodeclare Index mu, nu, rho, sig, chi, omega, xi, alpha, beta;
Autodeclare Index zz;

*-------------------------------------------------
* the momenta
Autodeclare Vector k,p,q;

*-------------------------------------------------
* A/B functions of sin/cos on the lattice
CF A, B, sf, AD, BD, del;
CF Subscript, T, Sin, Cos, SUM;

*-------------------------------------------------
* dirac matrices
F G, SIG;

*-------------------------------------------------
* integration constants and color functions
S FZERO, FONE, FTWO, FTHREE, rh, a;
CF ClrF(antisymmetric);


*---------------------------------------------------------------------------------
* GAUGE TENSOR FEYNMAN RULES, WITH SYMMETRIZATION, AND COLOR FACTORS


*==== FIRST ORDER FEYNMAN RULE ==*
L FirstOrderFeynmanRule = -rh/a*4*i_* 
                          (
                            FZERO  * T(0,1,0,p1,p2,mu) 
                          );

*==== SECOND ORDER FEYNMAN RULE ==*
L SecondOrderFeynmanRule = -rh/a*4*i_*
                           ( 
                            - FZERO * T(0,2,0,p1,p2,mu,nu)
                            + FONE  * ( T(1,0,1,p1,p2,mu,nu) + T(1,1,0,p1,p2,mu,nu) + T(0,1,1,p1,p2,mu,nu) )
                            - FTWO  * ( T(0,1,0,1,0,p1,p2,mu,nu) ) 
                           );

*==== THIRD ORDER FEYNMAN RULE ==*
L ThirdOrderFeynmanRule = -rh/a*4*i_*
                          (
                            - FONE * (  T(2,0,1,p1,p2,mu,nu,rho)  + T(2,1,0,p1,p2,mu,nu,rho) + T(1,0,2,p1,p2,mu,nu,rho)
                                       + T(1,2,0,p1,p2,mu,nu,rho) + T(1,1,1,p1,p2,mu,nu,rho) + T(0,2,1,p1,p2,mu,nu,rho)
                                       + T(0,1,2,p1,p2,mu,nu,rho) + T(0,3,0,p1,p2,mu,nu,rho) )

                             + FTWO * (  T(2,0,1,p1,p2,mu,nu,rho) + T(2,1,0,p1,p2,mu,nu,rho) + T(1,0,2,p1,p2,mu,nu,rho)
                                       + T(1,2,0,p1,p2,mu,nu,rho) + 2*T(1,1,1,p1,p2,mu,nu,rho) + T(0,2,1,p1,p2,mu,nu,rho)
                                       + T(0,1,2,p1,p2,mu,nu,rho)
                
                                       + T(1,0,1,0,1,p1,p2,mu,nu,rho) + T(1,0,1,1,0,p1,p2,mu,nu,rho) + T(1,1,0,1,0,p1,p2,mu,nu,rho)
	      		               + T(0,2,0,1,0,p1,p2,mu,nu,rho) + T(0,1,1,0,1,p1,p2,mu,nu,rho) + T(0,1,1,1,0,p1,p2,mu,nu,rho)
                                       + T(0,1,0,2,0,p1,p2,mu,nu,rho) + T(0,1,0,1,1,p1,p2,mu,nu,rho)
                                      )
			
			   + FTHREE * (
			                 T(2,0,1,p1,p2,mu,nu,rho) + T(2,1,0,p1,p2,mu,nu,rho) + T(1,2,0,p1,p2,mu,nu,rho)
                                       + T(1,1,1,p1,p2,mu,nu,rho) + T(0,2,1,p1,p2,mu,nu,rho) + T(0,1,2,p1,p2,mu,nu,rho)

				       + T(0,2,0,1,0,p1,p2,mu,nu,rho) + T(0,1,0,2,0,p1,p2,mu,nu,rho) 

                                       - T(0,1,0,1,0,1,0,p1,p2,mu,nu,rho)
                                      )
                                      
                          );



*===================================================
*  TRACE REPLACEMENTS
*===================================================

*set for the 1st order Feynman Rule.
*-----------------------------------------------------------------------------------------------------------------------
id T(0,1,0,p1?,p2?,mu?) = 1/(4*i_)*SIG(1,alpha,beta)
                                         *( A(0,p1,dum11)   *g_(1,dum1)  + B(0,p1,dum11) )
                                         *( AD(1,p1,p2,mu1) *g_(1,mu)    + BD(1,p1,p2,mu1) )
	                                 *( A(0,p2,dum22)   *g_(1,dum2)  + B(0,p2,dum22) );



************************************************************
* alternative fourier transform of the 2nd order 
* Feynman rule

*id T(0,2,0,p1?,p2?,mu?,nu?) = 1/(4*i_)*SIG(1,alpha,beta)
*                                         *( A(0,k1,dum11)  *g_(1,dum1)   + B(0,k1,dum11) )*
*                               del(mu,nu)*( AD(2,p1,p2,mu1)*g_(1,mu)     + BD(2,p1,p2,mu1))
*                                         *( A(0,k1,dum22)  *g_(1,dum2)   + B(0,k1,dum22) );


*id T(1,0,1,p1?,p2?,mu?,nu?) = 1/(4*i_)*SIG(1,alpha,beta)
*					*( A(1,p1,k1,mu1)  * g_(1,mu)    + B(1,p1,k1,mu1))
*					*( AD(0,k1,dum11)  * g_(1,dum1)  + BD(0,k1,dum11))
*					*( A(1,p2,k1,nu1)  * g_(1,nu)    + A(1,p2,k1,nu1));




*id T(1,1,0,p1?,p2?,mu?,nu?) = 1/(4*i_)*SIG(1,alpha,beta)
*                                         *( A(1,p1,k1,mu1)   *g_(1,mu)     + B(1,p1,k1,mu1) )
*					 *( AD(1,p2,k1,nu1)  *g_(1,nu)     + BD(1,p2,k1,nu1) )
*                                         *( A(0,k1,dum11)    *g_(1,dum1)   + B(0,k1,dum11) );


*id T(0,1,1,p1?,p2?,mu?,nu?) = 1/(4*i_)*SIG(1,alpha,beta)
*                                         *( A(0,k1,dum11)    *g_(1,dum1)   + B(0,k1,dum11) )
*                                         *( AD(1,p1,k1,mu1)  *g_(1,mu)     + BD(1,p1,k1,mu1) )
*                                         *( A(1,p2,k1,nu1)   *g_(1,nu)     + B(1,p2,k1,nu1) );



*id T(0,1,0,1,0,p1?,p2?,mu?,nu?) = 1/(4*i_)*SIG(1,alpha,beta)
*                                          *( A(0,k1,dum11)    *g_(1,dum1)   + B(0,k1,dum11) )
*                                          *( AD(1,p1,k1,mu1)  *g_(1,mu)     + BD(1,p1,k1,mu1) )
*                                          *( A(0,k1,dum22)    *g_(1,dum2)   + B(0,k1,dum22) )
*                                          *( AD(1,p2,k1,nu1)  *g_(1,nu)     + BD(1,p2,k1,nu1) )
*                                          *( A(0,k1,dum33)    *g_(1,dum3)   + B(0,k1,dum33) );



*************************************************************

* set for the 2nd order Feynman Rule.
*-----------------------------------------------------------------------------------------------------------------------
id T(0,2,0,p1?,p2?,mu?,nu?) = 1/(4*i_)*SIG(1,alpha,beta)
                                         *( A(0,p1,dum11)  *g_(1,dum1)   + B(0,p1,dum11) )*
                               del(mu,nu)*( AD(2,p1,p2,mu1)*g_(1,mu)     + BD(2,p1,p2,mu1))
                                         *( A(0,p2,dum22)  *g_(1,dum2)   + B(0,p2,dum22) );


id T(1,0,1,p1?,p2?,mu?,nu?) = 1/(4*i_)*SIG(1,alpha,beta)
                                         *( A(1,p1,k1,mu1)    *g_(1,mu)     + B(1,p1,k1,mu1) )
                                         *( AD(0,k1,dum11)    *g_(1,dum1)   + BD(0,k1,dum11) )
                                         *( A(1,k1,p2,nu1)    *g_(1,nu)     + B(1,k1,p2,nu1) );


id T(1,1,0,p1?,p2?,mu?,nu?) = 1/(4*i_)*SIG(1,alpha,beta)
                                         *( A(1,p1,k1,mu1)   *g_(1,mu)     + B(1,p1,k1,mu1) )
					 *( AD(1,k1,p2,nu1)  *g_(1,nu)     + BD(1,k1,p2,nu1) )
                                         *( A(0,p2,dum11)    *g_(1,dum1)   + B(0,p2,dum11) );


id T(0,1,1,p1?,p2?,mu?,nu?) = 1/(4*i_)*SIG(1,alpha,beta)
                                         *( A(0,p1,dum11)    *g_(1,dum1)   + B(0,p1,dum11) )
                                         *( AD(1,p1,k1,mu1)  *g_(1,mu)     + BD(1,p1,k1,mu1) )
                                         *( A(1,k1,p2,nu1)   *g_(1,nu)     + B(1,k1,p2,nu1) );


id T(0,1,0,1,0,p1?,p2?,mu?,nu?) = 1/(4*i_)*SIG(1,alpha,beta)
                                          *( A(0,p1,dum11)    *g_(1,dum1)   + B(0,p1,dum11) )
                                          *( AD(1,p1,k1,mu1)  *g_(1,mu)     + BD(1,p1,k1,mu1) )
                                          *( A(0,k1,dum22)    *g_(1,dum2)   + B(0,k1,dum22) )
                                          *( AD(1,k1,p2,nu1)  *g_(1,nu)     + BD(1,k1,p2,nu1) )
                                          *( A(0,p2,dum33)    *g_(1,dum3)   + B(0,p2,dum33) );


* set for the third order feynman rule
*-----------------------------------------------------------------------------------------------------------------------

id T(0,3,0,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
   				           *del(mu,nu)*del(mu,rho) *( A(0,p1,dum11)   * g_(1,dum1) + B(0,p1,dum11) )
                                                                   *( AD(3,p1,p2,mu1) * g_(1,mu)   + BD(3,p1,p2,mu1) )
                                                                   *( A(0,p2,dum22)   * g_(1,dum2) + B(0,p2,dum2));

id T(2,0,1,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
                                           *del(mu,nu)*( A(2,p1,k1,mu1) *g_(1,mu)     + B(2,p1,k1,mu1) )
                                                      *( AD(0,k1,dum11) *g_(1,dum1)   + BD(0,k1,dum11) )
                                                      *( A(1,k1,p2,rho1)*g_(1,rho)    + B(1,k1,p2,rho1) );

id T(2,1,0,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
                                           *del(mu,nu)*( A(2,p1,k1,mu1)  *g_(1,mu)    + B(2,p1,k1,mu1) )
                                                      *( AD(1,k1,p2,rho1)*g_(1,rho)   + BD(1,k1,p2,rho1) )
                                                      *( A(0,p2,dum11)  *g_(1,dum1)   + B(0,p2,dum11) );


id T(1,0,2,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
                                                      *( A(1,p1,k1,mu1) *g_(1,mu)    + B(1,p1,k1,mu1) )
                                                      *( AD(0,k1,dum11) *g_(1,dum1)  + BD(0,k1,dum1) )
                                          *del(nu,rho)*( A(2,k1,p2,rho1)*g_(1,rho)   + B(2,k1,p2,rho1) );


id T(1,2,0,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
                                                      *( A(1,p1,k1,mu1) *g_(1,mu)    + B(1,p1,k1,mu1) )
			                  *del(nu,rho)*( AD(2,k1,p2,rho1)*g_(1,rho)  + BD(2,k1,p2,rho1) )
                                                      *( A(0,p2,dum11) *g_(1,dum1)   + B(0,p2,dum1) );


id T(1,1,1,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
                                                      *( A(1,p1,k1,mu1) *g_(1,mu)    + B(1,p1,k1,mu1) )
                                                      *( AD(1,k1,k2,nu1)*g_(1,nu)    + BD(1,k1,k2,nu1) )
                                                      *( A(1,k2,p2,rho1) *g_(1,rho)  + B(1,k2,p2,rho1) );


id T(0,2,1,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
                                                      *( A(0,p1,dum11) *g_(1,dum1)   + B(0,p1,dum11) )
                                           *del(mu,nu)*( AD(2,p1,k1,mu1) *g_(1,mu)   + BD(2,p1,k1,mu1) )
                                                      *( A(1,k1,p2,rho1)*g_(1,rho)   + B(1,k2,p2,rho1) );


id T(0,1,2,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
                                                      *( A(0,p1,dum11) *g_(1,dum1)   + B(0,p1,dum11) )
                                                      *( AD(1,p1,k1,mu1) *g_(1,mu)   + BD(1,p1,k1,mu1) )
                                          *del(nu,rho)*( A(2,k1,p2,rho1)*g_(1,rho)   + B(2,k1,p2,rho1) );


id T(1,0,1,0,1,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
                                          *( A(1,p1,k1,mu1)   *g_(1,mu)     + B(1,p1,k1,mu1) )
                                          *( AD(0,k1,dum11)   *g_(1,dum1)   + BD(0,k1,dum11) )
                                          *( A(1,k1,k2,nu1)   *g_(1,nu)     + B(1,k1,k2,nu1) )
                                          *( AD(0,k2,dum22)   *g_(1,dum2)   + BD(0,k2,dum22) )
                                          *( A(1,k2,p2,rho1)  *g_(1,rho)    + B(1,k2,p2,rho1) );

id T(1,0,1,1,0,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
                                          *( A(1,p1,k1,mu1)   *g_(1,mu)     + B(1,p1,k1,mu1) )
                                          *( AD(0,k1,dum11)   *g_(1,dum1)   + BD(0,k1,dum1) )
                                          *( A(1,k1,k2,nu1)   *g_(1,nu)     + B(1,k1,k2,nu1) )
                                          *( AD(1,k2,p2,rho1) *g_(1,rho)    + BD(1,k2,p2,rho1) )
                                          *( A(0,p2,dum22)    *g_(1,dum2)   + B(0,p2,dum22) );

id T(1,1,0,1,0,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
                                          *( A(1,p1,k1,mu1)   *g_(1,mu)     + B(1,p1,k1,mu1) )
                                          *( AD(1,k1,k2,nu1)  *g_(1,nu)     + BD(1,k1,k2,nu1) )
                                          *( A(0,k2,dum11)    *g_(1,dum1)   + B(0,k2,dum11) )
                                          *( AD(1,k2,p2,rho1) *g_(1,rho)    + BD(1,k2,p2,rho1) )
                                          *( A(0,p2,dum22)    *g_(1,dum2)   + B(0,p2,dum22) );


id T(0,2,0,1,0,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
                                          *( A(0,p1,dum11)    *g_(1,dum1)   + B(0,p1,dum11) )
                               *del(mu,nu)*( AD(2,p1,k1,mu1)  *g_(1,mu)     + BD(1,p1,k1,mu1) )
                                          *( A(0,k1,dum22)    *g_(1,dum2)   + B(0,k1,dum22) )
                                          *( AD(1,k1,p2,rho1) *g_(1,rho)    + BD(1,k1,p2,rho1) )
                                          *( A(0,p2,dum33)    *g_(1,dum3)   + B(0,p2,dum33) );


id T(0,1,1,0,1,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
                                          *( A(0,p1,dum11)   *g_(1,dum1)   + B(0,p1,dum11) )
                                          *( AD(1,p1,k1,mu1) *g_(1,mu)     + BD(1,p1,k1,mu1) )
                                          *( A(1,k1,k2,nu1)  *g_(1,nu)     + B(1,k1,k2,nu1) )
                                          *( AD(0,k2,dum22)  *g_(1,dum2)   + BD(0,k2,dum22) )
                                          *( A(1,k2,p2,rho1) *g_(1,rho)    + B(1,k2,p2,rho1) );


id T(0,1,1,1,0,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
                                          *( A(0,p1,dum11)     *g_(1,dum1)   + B(0,p1,dum11) )
                                          *( AD(1,p1,k1,mu1)   *g_(1,mu)     + BD(1,p1,k1,mu1) )
                                          *( A(1,k1,k2,nu1)    *g_(1,nu)     + B(1,k1,k2,nu1) )
                                          *( AD(1,k2,p2,rho1) *g_(1,rho)   + BD(1,k2,p2,rho1) )
                                          *( A(0,p2,dum22)      *g_(1,dum2)    + B(0,p2,dum22) );


id T(0,1,0,2,0,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
                                          *( A(0,p1,dum11)   *g_(1,dum1)   + B(0,p1,dum11) )
                                          *( AD(1,p1,k1,mu1) *g_(1,mu)     + BD(1,p1,k1,mu1) )
                                          *( A(0,k1,dum22)   *g_(1,dum2)   + B(0,k1,dum22) )
                              *del(nu,rho)*( AD(2,k1,p2,nu1) *g_(1,nu)     + BD(2,k1,p2,nu1) )
                                          *( A(0,p2,dum33)   *g_(1,dum3)   + B(0,p2,dum33) );


id T(0,1,0,1,1,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
                                          *( A(0,p1,dum11)      *g_(1,dum1) + B(0,p1,dum11) )
                                          *( AD(1,p1,k1,mu1)    *g_(1,mu)   + BD(1,p1,k1,mu1) )
                                          *( A(0,k1,dum22)      *g_(1,dum2) + B(0,k1,dum22) )
                                          *( AD(1,k1,k2,nu1)    *g_(1,nu)   + BD(1,k1,k2,nu1) )
                                          *( A(1,k2,p2,rho1)    *g_(1,rho)  + B(1,k2,p2,rho1) );


id T(0,1,0,1,0,1,0,p1?,p2?,mu?,nu?,rho?) = 1/(4*i_)*SIG(1,alpha,beta)
                                          *( A(0,p1,dum11)      *g_(1,dum1) + B(0,p1,dum11) )
                                          *( AD(1,p1,k1,mu1)    *g_(1,mu)   + BD(1,p1,k1,mu1) )
                                          *( A(0,k1,dum22)      *g_(1,dum2) + B(0,k1,dum22) )
                                          *( AD(1,k1,k2,nu1)    *g_(1,nu)   + BD(1,k1,k2,nu1) )
                                          *( A(0,k2,dum33)      *g_(1,dum3) + B(0,k2,dum33) )
					  *( AD(1,k2,p2,rho1)   *g_(1,rho) + BD(1,k2,p2,rho1))
                                          *( A(0,p2,dum44)      *g_(1,dum4) + B(0,p2,dum44) );

*-----------------------------------------------------------------------------------------------------------------------
* == DEFINITION OF SIGMA ==
id SIG(1,zz1?,zz2?) = -i_/2 * ( g_(1,zz1,zz2) - g_(1,zz2,zz1) );


* == COMPUTE THE TRACE ==
trace4,1;

*------------------------------------------------------------------------------------------------------------------------

*id A(1,p1?,p2?,mu?)  = i_  * Cos(a * p1(mu)/2 + a * p2(mu));
*id B(1,p1?,p2?,mu?)  = rw  * Sin(a * p1(mu)/2 + a * p2(mu));
*id AD(1,p1?,p2?,mu?) = -i_ * Cos(a * p1(mu)/2 + a * p2(mu));
*id BD(1,p1?,p2?,mu?) = rw  * Sin(a * p1(mu)/2 + a * p2(mu));



* == REPLACE ALL A/B FACTORS ==
id A(0,p1?,mu?)  =   i_/a * SUM(mu) * Sin(a*p1(mu));
id B(0,p1?,mu?)  = 2*rw/a * SUM(mu) * Sin(a*p1(mu)/2)^2 - rh/a;
id AD(0,p1?,mu?) =  -i_/a * SUM(mu) * Sin(a*p1(mu));
id BD(0,p1?,mu?) = 2*rw/a * SUM(mu) * Sin(a*p1(mu)/2)^2 - rh/a;

id A(1,p1?,p2?,mu?)  = i_  * Cos(a * p1(mu)/2 + a * p2(mu)/2);
id B(1,p1?,p2?,mu?)  = rw  * Sin(a * p1(mu)/2 + a * p2(mu)/2);
id AD(1,p1?,p2?,mu?) = -i_ * Cos(a * p1(mu)/2 + a * p2(mu)/2);
id BD(1,p1?,p2?,mu?) = rw  * Sin(a * p1(mu)/2 + a * p2(mu)/2);

id A(2,p1?,p2?,mu?)   = -i_*a/2 * Sin(a * p1(mu)/2 + a * p2(mu)/2);
id B(2,p1?,p2?,mu?)   = a*rw/2  * Cos(a * p1(mu)/2 + a * p2(mu)/2);
id AD(2,p1?,p2?,mu?)  =  i_*a/2 * Sin(a * p1(mu)/2 + a * p2(mu)/2);
id BD(2,p1?,p2?,mu?)  = a*rw/2  * Cos(a * p1(mu)/2 + a * p2(mu)/2);

id A(3,p1?,p2?,mu?)   = -a^2/6  * Cos(a * p1(mu)/2 + a * p2(mu)/2);
id B(3,p1?,p2?,mu?)   = -a^2/6  * rw  * Sin(a * p1(mu)/2 + a * p2(mu)/2);
id AD(3,p1?,p2?,mu?)  = -a^2/6  * Cos(a * p1(mu)/2 + a * p2(mu)/2);
id BD(3,p1?,p2?,mu?)  = -a^2/6	* rw  * Sin(a * p1(mu)/2 + a * p2(mu)/2);


* == text parsing == *
*----------------------------------------
id d_(mu1?,mu2?) = del(mu1,mu2); 
.sort


*----------------------------------------
print FirstOrderFeynmanRule;

*print SecondOrderFeynmanRule;

*print ThirdOrderFeynmanRule;
*----------------------------------------

.end
