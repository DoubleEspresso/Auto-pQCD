(************************)
(*set the directory here*)
dir = SetDirectory["~/physics/research/perturbativeLattice/form_codes/v1/src/taylor.tmp"];

(*set the output dir to math.in*)

(******************************************)
(*  COMMON FEYNMAN RULE SUBSCRIPTING      *)
(* mathematica is surprisingly bad at     *)
(* appending  indices over several        *) 
(* momenta at once. e.g.                  *)
(*    Subscript[p+2q-3k,mu]               *)
(******************************************)

(*caveats :: cannot use momenta like p1,p2,p3*)
(*have to declare them differently*)
(*e.g. p,q,r,s,t etc.*)

subscript[input_, index_] := 
  Module[ {indexedAns, ans, expr, dir1, dir2},
 
    (*print to index.out in main directory tree*)  
    expr = ToString[input,FormatType->InputForm,CharacterEncoding->"ASCII"];
    Export["index.txt",expr];
  
    (*CALL PYTHON FROM HERE*)  
    indexedAns = ReadList[Evaluate["!python subscript.py index.txt -s " <> index]][[1]];
  Return[indexedAns];
];



(*a wrapper to ToString*)
tostring[mu_]:=ToString[mu, FormatType -> InputForm, CharacterEncoding -> "ASCII"];

(********************************************************************************)
(********************************************************************************)
(********************************************************************************)


(****************************************)
(*THE WILSON QCD VERTICES               *)
(****************************************)

(*first order qcd vertex*)
WQQG[f_,mu_]:=-(I*G[mu]*Cos[a*subscript[f,tostring[mu]]/2] + rw*Sin[a*subscript[f,tostring[mu]]/2]);

(*second order qcd vertex*)
WQQGG[f_,mu_,nu_]:=-a/2*d[mu,nu]*(-I*G[mu]*Sin[a*subscript[f,tostring[mu]]/2]+rw*Cos[a*subscript[f,tostring[mu]]/2]);




(*******************************************)
(* The Wilson Gauge action vertices        *)
(*******************************************)
WGGG[k1_, k2_, k3_, \[Mu]_, \[Nu]_, \[Lambda]_] := 
    Module[{ans, tmp1, tmp2, tmp3, idx1, idx2, idx3, i1Tab, i2Tab, i3Tab, pList, idxList},
			
     (*convert the indices to strings*)   
     idx1 = ToString[\[Mu], FormatType -> InputForm, CharacterEncoding -> "ASCII"];
     idx2 = ToString[\[Nu], FormatType -> InputForm, CharacterEncoding -> "ASCII"];
     idx3 = ToString[\[Lambda], FormatType -> InputForm, CharacterEncoding -> "ASCII"];
   
     tmp1 = subscript[Expand[k1 - k2], idx3];
     tmp2 = subscript[Expand[k2 - k3], idx1];
     tmp3 = subscript[Expand[k3 - k1], idx2];

  
    ans = -I*2/a*(d[\[Mu], \[Nu]]*Sin[(a tmp1)/2]*Cos[(a subscript[k3, idx1])/2] 
                + d[\[Nu], \[Lambda]]*Sin[(a tmp2)/2]*Cos[(a subscript[k1, idx2])/2] 
                + d[\[Lambda], \[Mu]]*Sin[(a tmp3)/2]*Cos[(a subscript[k2, idx3])/2]);
   
     Return[ans];
];


(******************************************)
(* THE WILSON PROPAGATORS                 *)
(******************************************)

(*quark propagator*)
WQQ[f_,mq_]:= a*(-I*SUM[tmp1]*G[tmp1]*Sin[a*subscript[f,"tmp1"]] + a*mq + 2*rw*SUM[tmp1]*Sin[a*subscript[f,"tmp1"]/2]^2)/
  DQ[ SUM[tmp2]*Sin[a*subscript[f,"tmp2"]]^2 + 4*rw^2*SUM[tmp2]*Sin[a*subscript[f,"tmp2"]/2]*SUM[tmp3]*Sin[a*subscript[f,"tmp3"]/2]^2+4*rw*a*mq*SUM[tmp2]*Sin[a*subscript[f,"tmp2"]/2]^2 + a^2*mq^2 ];

(*gluon propagator*)
(*note that the gauge parameter has been set to 1 by default*)
WGG[f_,mu_,nu_]:= a^2 * d[mu,nu]/DB[4*SUM[tmp1]*Sin[a*subscript[f,"tmp1"]/2]^2];


(********************************************************************************)
(********************************************************************************)
(********************************)
(*  THE OVERLAP FENMAN RULES    *)
(********************************)




(********************************************************************************)
(********************************************************************************)
(*************************************)
(*  THE QUARK ANGULAR MOMENTUM OPS   *)
(*  color is computed by hand here   *)
(*************************************)
QOP1[f_]:=I/a*SUM[tmp1]*SUM[tmp2]*q[tmp1]*q[tmp2]*G[tmp1]*Sin[a*subscript[f,"tmp2"]];

(*f is the sum of incoming and outgoing quark momenta for this graph*)
QOP2[f_]:=I*SUM[tmp1]*SUM[tmp2]*q[tmp1]*q[tmp2]*G[tmp1]*Cos[a*subscript[f,"tmp2"]/2];




(******************************************)
(*  TAYLOR EXPANSION ROUTINES             *)
(*  Currently the maximum order is 2      *)
(*  significant changes must be made to   *)
(*  extend this to higher orders          *)
(******************************************)


(************************************************)
(* Non-commutative expansion is done using FORM *)
(* these codes are necessary to properly taylor *)
(* the overlap vertices                         *)
(************************************************)
NCExpand[ansExpr_]:= Module[ {FormOutput,ans,expr,dir1,dir2},
    

  expr=ToString[ansExpr,FormatType->InputForm,CharacterEncoding->"ASCII"];
  Export[dir<>"/ncmExpand.txt",expr];

  (*CALL PYTHON FROM HERE*)
  FormOutput=ReadList[Evaluate["!sh lattice_expand.sh"]][[1]];

  (*replacements*)
  ans=FormOutput//.{k[mu_]->Subscript[k, mu],p[mu_]->Subscript[p, mu],g[x__]->G[x]};

  Return[ans];
]



(*************************************)
(* The actual taylor expansion codes *)
(* these work up to 2nd order in the *)
(* the external momentum             *)
(*************************************)


(*FIRST ORDER TAYLOR EXPANSION*)
LatticeTaylorExpand1[expr_,p_,idxList_]:=
  Module[
    
   {zeroD,firstD,secondD,repZero,
    tmp1,tmp2,tmp3,tmp4,tmp5,trigReps,
    qpropRep, gpropRep, ans2, ans3},


    zeroD  = Expand[expr];
    firstD = D[Expand[expr],p]//.{Derivative[1, 0][Subscript][p, nu_]->d[nu,idxList]};

    repZero = {Subscript[p,nu_]->0,p->0,Subscript[0,nu_]->0};

    
    trigReps = {
      Sin[a*(-Subscript[k,mu_])]->-B[mu],
      Cos[1/2*a*(-Subscript[k,mu_])]->M[mu],
      Sin[1/2*a*(-Subscript[k,nu_])]->-S[nu],
      Cos[a*Subscript[k,nu_]/2]->M[nu],
      Sin[1/2*a*Subscript[k,nu_]]->S[nu],
      Sin[a*Subscript[k,nu_]]->B[nu],
      Cos[a*Subscript[k,nu_]]->M[nu]^2-S[nu]^2,
      SUM[nu_]^2 -> SUM[nu]
    };


    qpropRep = {Derivative[1][DQ][x_] -> 1, DQ[x__] -> DQ};
    gpropRep = {Derivative[1][DB][x_] -> 1, DB[x__] -> DB};

    tmp1 = firstD//.repZero;
    tmp2 = zeroD//.repZero;

    tmp3 = tmp1//.trigReps;
    tmp4 = tmp2//.trigReps;
    
    ans = tmp4 + SUM[idxList]*p[idxList]*tmp3;

   ans2 = ans//.qpropRep;
   ans3 = ans2//.gpropRep;

  Return[ans3];
]

(****************************************)
(* SECOND ORDER TAYLOR EXPANSION        *)
(* just one external  momentum          *)
(****************************************)
LatticeTaylorExpand2[expr_,p_,sum1_,sum2_]:=
  Module[
	 {DERIV1,DERIV2,DERIVREP0,ZERODERIVREP0,
	  DERIV2REP0,FINALREP,FINALREP2,FINALREP0,
	  ZERODERIV,ZERODERIV2,ZERODERIV3,PIECEZERO,
	  PIECEONE,PIECETWO,ans},

	 (*COMPUTE THE FIRST DERIVATIVE*)
	 DERIV1=D[Expand[expr],p]//.{ Derivative[1, 0][Subscript][p, nu_]->d[nu,sum1],Derivative[1, 0][Subscript][-p, nu_]->-d[nu,sum1]};
	 DERIV2=D[Expand[DERIV1],p]//.{ Derivative[1, 0][Subscript][p, nu_]->d[nu,sum2],Derivative[1, 0][Subscript][-p, nu_]->-d[nu,sum2]};

	 (*add to the first derivative, the piece with p->0*)
	 ZERODERIV=expr;

	 (*MAKE REPLACEMENTS IN ZERO DERIVATIVE*)
	 FINALREP0={Subscript[p, nu_]->0,p->0,Subscript[0, nu_]->0};
	 ZERODERIVREP0=ZERODERIV//.FINALREP0;

	 (*MAKE REPLACEMENTS IN FIRST DERIVATIVE*)
	 FINALREP0={Subscript[p, nu_]->0,p->0,Subscript[0, nu_]->0};
	 DERIVREP0=DERIV1//.FINALREP0;

	 (*MAKE REPLACEMENTS IN SECOND DERIVATIVE*)
	 FINALREP0={Subscript[p, nu_]->0,p->0,Subscript[0, nu_]->0};
	 DERIV2REP0=DERIV2//.FINALREP0;

         FINALREP2 = {
           Sin[a*(-Subscript[k,mu_])]->-B[mu],
           Cos[1/2*a*(-Subscript[k,mu_])]->M[mu],
           Sin[1/2*a*(-Subscript[k,nu_])]->-S[nu],
           Cos[a*Subscript[k,nu_]/2]->M[nu],
           Sin[1/2*a*Subscript[k,nu_]]->S[nu],
           Sin[a*Subscript[k,nu_]]->B[nu],
           Cos[a*Subscript[k,nu_]]->M[nu]^2-S[nu]^2,
           SUM[nu_]^2 -> SUM[nu]
         };



	 derivReps = {Derivative[1][DQ][x_] -> 1, Derivative[2][DQ][x_]->0, DQ[x__] -> DQ,
	              Derivative[1][DB][x_] -> 1, Derivative[2][DB][x_]->0, DB[x__] -> DB};


	 PIECEONE=DERIVREP0//.FINALREP2;
	 PIECEZERO=ZERODERIVREP0//.FINALREP2;
	 PIECETWO=DERIV2REP0//.FINALREP2;

	 ans = PIECEZERO+SUM[sum1]*p[sum1]*PIECEONE+1/2*SUM[sum1]*SUM[sum2]*p[sum1]*p[sum2]*PIECETWO;

	 Return[ans//.derivReps];
]

(*****************************************************)
(*****************************************************)
(*****************************************************)





dir=SetDirectory["~/physics/research/perturbativeLattice/form_codes/v1/src/taylor.tmp"]
printDir="~/physics/research/perturbativeLattice/form_codes/v1/src/";

(*printing script*)
PrintToFile[expr_]:= Module[{ans,file}, 
  ans=ToString[expr, FormatType -> InputForm, CharacterEncoding -> "ASCII"];
  Export[printDir<>"tst.txt", ans]];

 v0 = LatticeTaylorExpand2[WQQ[k+p,0]/.{tmp1->ix3,tmp2->ix4,tmp3->ix5},p,ix1,ix2]//Expand;
v1 = LatticeTaylorExpand2[WQQ[k,0]/.{tmp1->ix8,tmp2->ix9,tmp3->ix10},p,ix6,ix7]//Expand;

expr=1*1/a^4*v0.v1;

PrintToFile[expr];
Exit[]
