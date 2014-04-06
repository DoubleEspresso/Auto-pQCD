#procedure expandSumsOverlap()

* // this routine will expand the sums inserted in the monomial
* // in python, and prep the expression for the IR-reduce procedure 
* // in FORM.
* //=================================================================

* // we first re-write the powers of kh in the expression
* // and format them into F-constants
repeat;
  id kh(1) = x1;
  id kh(2) = x2;
  id kh(3) = x3;
  id kh(4) = x4;
endrepeat;

* // make the f-function constants
id x1^pow? = f(pow);
id x2^pow? = f(pow);
id x3^pow? = f(pow);
id x4^pow? = f(pow);



id f(x1?)*f(x2?)*f(x3?)*f(x4?) = FS(x1/2,x2/2,x3/2,x4/2);
id f(x1?)*f(x2?)*f(x3?)        = FS(x1/2,x2/2,x3/2,0);
id f(x1?)*f(x2?)               = FS(x1/2,x2/2,0,0);
id f(x1?)                      = FS(x1/2,0,0,0);

* //========================================================================
* // include the various powers of quark/gluon propagators into the 
* // integration shorthands.  Note that the power on the boson propagator 
* // is listed first.  We also multiply them by (-1) to absorb the sign.
id FS(pw1?,pw2?,pw3?,pw4?)*db^pw5?*dq^pw6?*aa^pw7?*bb^pw8? = F(-pw6,-pw5,pw1,pw2,pw3,pw4,pw7,pw8);
