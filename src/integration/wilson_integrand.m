#!/usr/local/bin/MathematicaScript -script

(***************************************)
(* read in usr args*)
(***************************************)
(*err=ToExpression[$ScriptCommandLine[[2]]];*)
err = 9;

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


(*THE INTEGRAND*)
integrand = If[fOrj==0, kh1^nt*kh2^nx*kh3^ny*kh4^nz/(dw^powq dg^powg), 
	       kh1^nt*kh2^nx*kh3^ny*kh4^nz/(dw^powq dg^powg) - kh1^nt*kh2^nx*kh3^ny*kh4^nz/( dg^(powq+powg))];


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
pows=StringJoin[ToString[powq],",",ToString[powg],",",ToString[nx],",",ToString[ny],",",ToString[nz],",",ToString[nt],"]->"]
intStr=StringJoin[head,pows,ToString[SetPrecision[ans,err],FormatType -> InputForm, CharacterEncoding ->"ASCII" ],","];


(*it is actually difficult to get python to read results from mathematica*)
(*so we export the result to a file, and have the shell print out the result*)
intStr >>> math.txt;

Exit[]; 

