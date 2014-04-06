integrand = pf*J[1,2,0,0,0,1];
rw = 1.0;
rho= 1.0;

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

(*set rho here*)
rho = 1;

dov = 2 * rho * Sqrt[sumk + (2 * rw * sumkhalf - rho)^2] + 2 * rho * (2 * rw * sumkhalf - rho);
wk = Sqrt[sumk + (2 * rw * sumkhalf - rho)^2];
aa = 1/wk;
bb = 1/(wk+rho);


(*THE INTEGRAND*)
J[q_,b_,nt_,nx_,ny_,nz_]:=kh1^nt*kh2^nx*kh3^ny*kh4^nz/(dw^q dg^b) - kh1^nt*kh2^nx*kh3^ny*kh4^nz/( dg^(b+q))

F[q_,b_,nt_,nx_,ny_,nz_]:=kh1^nt*kh2^nx*kh3^ny*kh4^nz/(dw^q dg^b) 

(*INTEGRATION ROUTINE*)
ans = NIntegrate[integrand,
		 {k1,-Pi,Pi},
		 {k2,-Pi,Pi},
		 {k3,-Pi,Pi},
		 {k4,-Pi,Pi},
		 MaxRecursion->30,
		 Method->{"LocalAdaptive",Method->"ClenshawCurtisOscillatoryRule"},
		 PrecisionGoal->10];


Print[ToString[SetPrecision[ans,11], FormatType -> InputForm, 
	 CharacterEncoding -> "ASCII"]];

Exit[]; 