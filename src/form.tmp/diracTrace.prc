#procedure diracTrace()


* // traces of up to six gamma matrices are computed by this routine.
* // These replacements have been pre-computed using FeynCalc in D-dimensions.
*==========================================================================


*// case when no dirac matrices are present
*// we introduce a scale factor to catch this
id G(?dum1) = sf * G(?dum1);

* // ======================================================================
* //
* //  make reps for gamma matrices
* //

id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?)*G(dum5?)*G(dum6?)*G(dum7?)*G(dum8?)*G(dum9?)*G(dum10?) = ORDERTEN;
id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?)*G(dum5?)*G(dum6?)*G(dum7?)*G(dum8?)*G(dum9?) = 0;
id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?)*G(dum5?)*G(dum6?)*G(dum7?)*G(dum8?) = ORDEREIGHT;
id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?)*G(dum5?)*G(dum6?)*G(dum7?) = 0;
id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?)*G(dum5?)*G(dum6?) = G(dum1,dum2,dum3,dum4,dum5,dum6);
id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?)*G(dum5?) = 0;
id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?) = G(dum1,dum2,dum3,dum4);
id G(dum1?)*G(dum2?)*G(dum3?) = 0;
id G(dum1?)*G(dum2?) = G(dum1,dum2);


  
*// ORDER SIX REPLACEMENTS
id G(dum1?,dum2?,dum3?,dum4?,dum5?,dum6?) =     
    4*ld(dum1, dum6)*ld(dum2, dum5)*ld(dum3, dum4) - 
    4*ld(dum1, dum5)*ld(dum2, dum6)*ld(dum3, dum4) - 
    4*ld(dum1, dum6)*ld(dum2, dum4)*ld(dum3, dum5) + 
    4*ld(dum1, dum4)*ld(dum2, dum6)*ld(dum3, dum5) + 
    4*ld(dum1, dum5)*ld(dum2, dum4)*ld(dum3, dum6) - 
    4*ld(dum1, dum4)*ld(dum2, dum5)*ld(dum3, dum6) + 
    4*ld(dum1, dum6)*ld(dum2, dum3)*ld(dum4, dum5) - 
    4*ld(dum1, dum3)*ld(dum2, dum6)*ld(dum4, dum5) + 
    4*ld(dum1, dum2)*ld(dum3, dum6)*ld(dum4, dum5) - 
    4*ld(dum1, dum5)*ld(dum2, dum3)*ld(dum4, dum6) + 
    4*ld(dum1, dum3)*ld(dum2, dum5)*ld(dum4, dum6) - 
    4*ld(dum1, dum2)*ld(dum3, dum5)*ld(dum4, dum6) + 
    4*ld(dum1, dum4)*ld(dum2, dum3)*ld(dum5, dum6) - 
    4*ld(dum1, dum3)*ld(dum2, dum4)*ld(dum5, dum6) + 
    4*ld(dum1, dum2)*ld(dum3, dum4)*ld(dum5, dum6);


*// ORDER FOUR REPLACEMENTS
id G(dum1?,dum2?,dum3?,dum4?) = 
   4*(ld(dum1, dum4)*ld(dum2, dum3) - 
   ld(dum1, dum3)*ld(dum2, dum4) + 
   ld(dum1, dum2)*ld(dum3, dum4));


*// ORDER TWO REPLACEMENTS
id G(dum1?,dum2?) = 4*ld(dum1, dum2);



*// terms second order in sf contained dirac matrices
*// but terms linear in sf did not, they are replaced
*// with the D (dimension) --> 4

id sf^2 = 1;
id sf = 4;



*//
*// after these replacements, we make gamma a commuting function
*// so that the contract delta routines work properly
*//
id G(dum1?) = g(dum1);

*================================================================
* for un-polarized renormalization calculations, twist-2 operators
* etc. the sigma-tensor will vanish always, we make these replacements
* here.

*------------------------------
*/ reps for sig-tensor
*/ ensure that sig is anti-symmetric
*--------------------------------
id sig(dum1?,dum1?) = 0;
id sig(dum1?,dum2?)*p1?(dum1?)*p1?(dum2?) = 0;
id sig(dum1?,dum2?)*m?(dum1?)*m?(dum2?) = 0;
id sig(dum1?,dum2?)*s?(dum1?)*s?(dum2?) = 0;
id sig(dum1?,dum2?) = 0;

*--------------------
* FOR NOW SET EPS TO ZERO
*---------------------
id eps(dum1?,dum2?,dum3?,dum4?) = 0;


*---------------------------
* RESET THE ld(dum1?,dum2?)
*---------------------------
repeat;
  id ld(dum1?,dum2?) = d(dum1,dum2);
endrepeat;
