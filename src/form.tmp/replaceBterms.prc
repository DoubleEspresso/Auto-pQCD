#procedure replaceBterms()

* =========================================
* INTEGRATION ROUTINE
* =========================================
* //
* // remove the various b-factors in terms of
* // their definitions
* //
id b(dum1?)*b(dum1?) = 4 * s(dum1)^2 - 4*s(dum1)^4;


* //
* // set all odd powers of b to zero
* // a 2nd time just to be sure..
* //
id b(dum1?) = sf*b(dum1);
id sf^pow?!even_ =0;
id sf = 1;

* //
* // next replace powers of s(x) --> 1/2 * kh(x)
* //
repeat;
  id s(dum1?) = 1/2*kh(dum1);
endrepeat;


* // at this point we should have replaced all B-terms with their corresponding
* // kh-equivalents.  If any b-terms are left, something went wrong,
id b(dum1?) = ErrorFoundBTermPresent(dum1);
id ErrorFoundBTermPresent(dum1?) = 0;