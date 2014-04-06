#!/usr/local/bin/MathematicaScript -script

(***************************************)
(* read in usr args*)
(***************************************)
(*err=ToExpression[$ScriptCommandLine[[2]]];*)
err = 6;

powq=ToExpression[$ScriptCommandLine[[3]]];
powg=ToExpression[$ScriptCommandLine[[4]]];

nx=ToExpression[$ScriptCommandLine[[5]]];
ny=ToExpression[$ScriptCommandLine[[6]]];
nz=ToExpression[$ScriptCommandLine[[7]]];
nt=ToExpression[$ScriptCommandLine[[8]]];

rw=ToExpression[$ScriptCommandLine[[9]]];
rho=ToExpression[$ScriptCommandLine[[10]]];
fOrj=ToExpression[$ScriptCommandLine[[11]]];

coreNum=ToExpression[$ScriptCommandLine[[12]]];

aPow = ToExpression[$ScriptCommandLine[[13]]];
bPow = ToExpression[$ScriptCommandLine[[14]]];

Clear[Coeff];

(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
(*Common definitions for denominators*)
(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
sumk = Sin[k1]^2 + Sin[k2]^2 + Sin[k3]^2 + Sin[k4]^2;
sumkhalf = Sin[k1/2]^2 + Sin[k2/2]^2 + Sin[k3/2]^2 + Sin[k4/2]^2;

(*Commond definitions for khat variables*)
kh1 = 4*Sin[k1/2]^2;
kh2 = 4*Sin[k2/2]^2;
kh3 = 4*Sin[k3/2]^2;
kh4 = 4*Sin[k4/2]^2;

(*prefactor for all lattice integrands*)
pf = 1/(2*Pi)^4;


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
(*DEFINITIONS FOR WILSON FERMIONS*)
(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
dg = 4 * sumkhalf;
dw = sumk + (2 * rw * sumkhalf)^2;


(*%%%%%%%%%%%%%%%%%%%%%%*)
(*DEFINITIONS FOR OVERLAP FERMIONS*)
(*%%%%%%%%%%%%%%%%%%%%%%%*)
dov = 2 * rho * Sqrt[sumk + (2 * rw * sumkhalf - rho)^2] + 2 * rho * (2 * rw * sumkhalf - rho);
wk = Sqrt[sumk + (2 * rw * sumkhalf - rho)^2];
aa = 1/wk;
bb = 1/(wk+rho);

(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
(* NOTE: for IR divergent integrals          *)
(* replace one fermion prop with             *)
(* 1/q = 1/g + (1/q-1/g)                     *)
(* reduces the degree of divergence by two   *)
(* for wilson, by one for overlap.           *)
(*                                           *)
(* integrals are IR divergent by a simple    *)
(* counting                                  *)
(*                                           *)
(* INT = k_1^{2n1} k_2^{2n2} k_3^{2n3}       *) 
(*       k_4^{2n4} / (dq^nq * dg^ng)         *)
(* if nq + ng < 2 + n1 + n2 + n3 + n4        *)
(*      --> IR SAFE                          *)
(*                                            *)
(* for overlap there are a couple differences:*)
(* there are two additional factors to        *)
(* consider ::                                *)
(*    aa = wk^-1 and bb = (wk + rho)^-1       *) 
(* luckily powers of these do not affect      *)
(* the power counting, which                  *)
(* follows wilson fermions                    *)
(* with the exception that the reduction      *)
(* procedure reduces the IR divergence by 1   *)
(* instead of 2 and that powers of bb are     *)
(* accompanied by powers of 1/2 e.g.          *)
(* 1/(dq*dg)*bb**2 - 1/dg**2 * 1/2**2 etc.    *)
(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
(* Script to compute Coeff -- making sure the  *)
(* integration is convergent in the k->0 limit *)
(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
SolveCoeff[integrand_, powList_] := 
 Module[{expr, tmp, tmp2, tmp3, kToZero, kToExpand, kExpand, kLimit, 
   zeroCount, nonZeroCount, reps, nReps, c},
  (*initialize*)
  zeroCount = 0;
  nonZeroCount = 0;
  nReps = "";
  
  (*collect the zero/nonzero k's appearing in the numerator*)

  Do[
   If[powList[[i]] == 0,
    zeroCount++;
    kToZero[zeroCount] = i;
    ];
   
   If[powList[[i]] != 0,
    nonZeroCount++;
    kToExpand[nonZeroCount] = i;
    ];
   , {i, 4}];
  
  (*case when all powers of k are present*)
  If[zeroCount == 0,
   kExpand[1] = k1; kExpand[2] = k2; kExpand[3] = k3; kLimit = k4;
   nReps = "{}"];
  
  (*case when no powers of k are present*)
  If[zeroCount == 4,
   kLimit = k1;
   nReps = ToString[{k2 -> 0, k3 -> 0, k4 -> 0}]];
  
  (*case when some powers of k are present*)
  Do[
   If[zeroCount != 0 && zeroCount < 4, 
    nReps = StringJoin[nReps, "k", ToString[kToZero[i]], "->0,"]], {i,
     zeroCount}
   ];
  
  (*the taylor expanded momenta*)
  Do[
   If[nonZeroCount > 0 && nonZeroCount < 4, 
    kExpand[i] = 
     ToExpression[StringJoin["k", ToString[kToExpand[i]]]]], {i, 
    nonZeroCount - 1}
   ];
  

  (*the momenta limit momenta->0*)  
  If[nonZeroCount > 0 && nonZeroCount < 4, 
   kLimit = 
    ToExpression[
     StringJoin["k", ToString[kToExpand[nonZeroCount]]]]];
  
  (*make the replacement rules*)
  
  If[zeroCount != 0 && zeroCount < 4, 
   nReps = StringJoin["{", StringDrop[nReps, -1], "}"]];
  nReps = ToExpression[nReps];
  
  
  (*debug*)
  Print[nReps];
  Print["Expansion Variables"];
  Do[Print[kExpand[i]],{i,nonZeroCount-1}];
  Print["Limit Variables"];
  Print[kLimit];




  partCoeff    = Coefficient[integrand,Coeff];
  partNonCoeff = Coefficient[integrand*Coeff,Coeff];

 Print["the integrand", integrand];
 Print["the coeff part",partCoeff];
 Print["the non-coeff part", partNonCoeff];

 (*  Print["coefficient term, ",partCoeff*Coeff];
  Print["non-coefficient term, ", partNonCoeff];*)





  (* avoid calling solve directly*)
  byHand = (-partNonCoeff/partCoeff)/.nReps;


  Print[byHand];
  Print[kLimit];




  (*byHand = Limit[byHand,kLimit->0];

  Print[byHand];

  tmp = integrand /. nReps;

  factor = 1; (*(dw^powq dg^powg)*(dg^(powq+powg))/.nReps;*)


  Print["Solving for Coeff"];
  Print[tmp];


  (*tmp2 = Solve[(dw^powq dg^powg) (dg^(powq + powg)) tmp == 0, 
     Coeff][[1, 1]];*)

  tmp2 = Solve[factor * tmp == 0, 
     Coeff][[1, 1]];

  tmp3 = Coeff /. tmp2;*)







  (*taylor expansion*)
  
  (*Do[tmp3 = Normal[Series[tmp3, {kExpand[i], 0, 1}]], {i, 
    nonZeroCount - 1}];*)


 If[nonZeroCount > 1, Do[byHand = Normal[Series[byHand, {kExpand[i], 0, 1}]], {i, 
    nonZeroCount - 1}]];
  

  (*limit, was tmp3*)
  c = Limit[byHand, kLimit -> 0];

  
  Return[c];
  ]



(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
(*THE INTEGRAND                    *)
(* this form of the integrand is   *)
(* used for glue angular momentum  *)
(* calculations.  We do not use dov*)
(* in this code!!                  *)
(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

(* NB : for quark/glue angular momentum operators, we use dw for the quark propagator!!*)

integrand = If[fOrj==0, kh1^nt*kh2^nx*kh3^ny*kh4^nz/(dw^powq dg^powg)*aa^aPow * bb^bPow, 
	       kh1^nt*kh2^nx*kh3^ny*kh4^nz/(dw^powq dg^powg + eps)*aa^aPow * bb^bPow - Coeff * kh1^nt*kh2^nx*kh3^ny*kh4^nz/( dg^(powq+powg) + eps)];
 

(* The coefficient must now be fixed to ensure  *)
(* that the integral converges for k->0         *)
(* to accomplish this, we taylor-expand         *) 
(* the J-integrand around k=0 and solve         *)
(* numerically the condition that               *)
(*     i1 - Coeff*i2 = 0,                       *)
(* which ensures a finite result for the        *)
(* integrand. We return this coeff to the       *)
(* main program so it can be used for the       *)
(* bosonic reduction calculation as well.       *)

(*Coeff = If[fOrj==0, 1, SolveCoeff[integrand,{nt,nx,ny,nz}]];*)
Coeff = If[fOrj==0, 1, 0.25];
eps = 1*10^(-10);


(*INTEGRATION ROUTINE*)
ans = NIntegrate[pf * integrand,
		 {k1,-Pi,Pi},
		 {k2,-Pi,Pi},
		 {k3,-Pi,Pi},
		 {k4,-Pi,Pi},
		 MaxRecursion->30,
		 Method->{"LocalAdaptive",Method->"ClenshawCurtisOscillatoryRule"},
		 PrecisionGoal->err];



(*set a string*)
head=If[fOrj==0,"F[","J["];
pows=StringJoin[ToString[powq],",",ToString[powg],",",ToString[nx],",",ToString[ny],",",ToString[nz],",",ToString[nt],",",ToString[aPow],",",ToString[bPow],"]->"];
intStr=StringJoin[head,pows,ToString[SetPrecision[ans,err],FormatType -> InputForm, CharacterEncoding ->"ASCII" ],","];

cArg = StringJoin[ToString[powq],",",ToString[powg],",",ToString[nx],",",ToString[ny],",",ToString[nz],",",ToString[nt],",",ToString[aPow],",",ToString[bPow]];

coeffStr=If[fOrj==0,"1",StringJoin["c[",cArg,"]->",ToString[Coeff],","]];

If[fOrj==1,Print[coeffStr]];

(*it is actually difficult to get python to read results from mathematica*)
(*so we export the result to a file, and have the shell print out the result*)
intStr >>> math.txt;
If[fOrj==1,coeffStr >>> math.txt];

Exit[]; 
