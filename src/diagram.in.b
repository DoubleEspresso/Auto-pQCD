# this file controls the input python and FORM will see
# it specifies options, variables and calculation techniques.
#----------------------------------------------------------------------------------------------------------------------------------
# variables include:
#   o the explicit value any external lorentz indices.
#   o the order of the Taylor expansion.
#   o the final operator structure to be projected out.
#   o the method used to numerically integrate (mathematica or VEGAS).
#   o whether traces should be performed over all dirac matrices
#   o the diagram to be computed.
#
#
#   to specify a diagram, we specify the vertices, and propagators, along with lorentz indices
#   and momentum flow.  There is little error checking in the codes, so be careful.  The overall
#   factors of 1/a^4 need to be included by the user as well.  The formatting for various vertices and
#   propagators is specified below ::
#
#-----------------------------------------------------------------------------------------------------------------------------------
#   Wilson Action:
#
#   QCD - Vertices :: only up to order g0^2 is implemented, and only 2nd order Taylor expansions.
#                     for all vertices, plus momenta is flowing inward!  Color computations are not 
#                     implemented yet!!
#
#       -- 1pt Vertex --          WQQG[f,index],  where p2 = p1 + k1, k1 being the gluon momentum flowing into the vertex.
#
#       -- 2pt Vertex --          WQQGG[f,index1,index2].      
#
#
#       -GLUON VERTICES-
#
#       -- 3pt Vertex --         WGGG[p1,p2,p3,index1,index2,index3,c1,c2,c3], (+momentum assigned in and CW)
#
#       -- 4pt Vertex --         WGGGG[p1,p2,p3,p4,i1,i2,i3,i4], (this is actually split into two pieces)
#
#
#
#       -PROPAGATORS-
#       -- Quark --              WQQ[p1,m0]
#       -- Gluon --              WGG[p1,ix1,ix2]
#
#------------------------------------------------------------------------------------------------------------------------
#
#  Quark/Gluon angular momentum operators
#  
#  --Quark Operator-- (twist-2,spin-2,contracted with light-like vector to project out the symmetrized and traceless piece)
#
#    order g0:    QOP1[p1,p2,ex1,ex2]  p1 = incoming quark, p2 = outgoing quark, ex1,ex2 are external indices.
#
#    order g0^2:  QOP2[p1,p2,ex1,ex2]
#
#
#  --Gluon Operator-- (twist-2, spin-2 defined from the overlap derivative)
#
#    order g0:   GOP1[p1,p2,ex1,ex2] 
#
#    order g0^2  GOP2[p1,p2,p3,ex1,ex2,ex3]
#
#    order g0^3  GOP3[p1,p2,p3,p4,ex1,ex2,ex3,ex4], color computed separately.
#
#--------------------------------------------------------------------------------------------------------------------------
#
#
#  Example input file:
#
#  $$ OutFile     :: QuarkWFR.Wilson
#  $$ TaylorOrder :: 1
#  $$ Method      :: math
#  $$ Trace       :: false
#  $$ External    :: none
#  $$ Sum         :: ix1,ix2
#  $$ Factor      :: 1/a^4
#  $$ Operator    :: (i_)*sum(dum1?)*p1(dum1?)*g(dum1?) // this should be formatted for FORM.
#  $$ Diagram     :: WQQG[p,k,ix1].WQQ[k,0].WQQG[k,p,ix2].WGG[k]
#
#
#
#  pre-checked diagrams ::
#  
#  for quark-wfr (Wilson action):
#  rainbow diagram :: $$ Diagram       ::  WQQG[p+k,ix1].WQQ[k,0].WQQG[k+p,ix2].WGG[k-p,ix1,ix2]
#  tadpole diagram :: 
#
#  for twist-2 mixing (Q->Q):
#  vertex-diagram ::  $$ Diagram       ::  WQQG[p+k,ix1].WQQ[k,0].QOP1[k,ex1,ex2].WQQ[k,0].WQQG[k+p,ix2].WGG[k-p,ix1,ix2]
#  sails-diagrams ::  $$ Diagram       ::  
#
#  for twist-2 mixing (Q->G): (2nd order T. expansion required, with dirac-trace)
#  quark-loop vertex :: $$ Diagram     ::  WQQG[2*k+p,ex1].WQQ[k,0].QOP1[k,ex3,ex4].WQQ[k,0].WQQG[2*k+p,ex2].WQQ[k+p,0]
#
#
#  this will compute the quark wave-function renormalization diagram, and collect only those terms which have 
#  momentum structure given by the Operator line.  Printing all results to QuarkWFR.Wilson.
#----------------------------------------------------------------------------------------------------------------------------
#
# $$ BEGIN DIAGRAM.IN

$$ OutFile       ::  wilson.quark.wfr.txt
$$ Action        ::  Overlap
$$ TaylorOrder   ::  1
$$ Method        ::  mathematica
$$ Integrate     ::  true
$$ Error         ::  6
$$ Trace         ::  false
$$ rw            ::  0.1
$$ rho           ::  1.0
$$ Sum           ::  
$$ ExtIndices    ::  
$$ ExtMomentum   ::  p
$$ Factor        ::  1/a^4
$$ Diagram       ::
#$$ Diagram       ::  WQQG[2*k+p,ex1].WQQ[k,0].QOP1[k].WQQ[k,0].WQQG[2*k+p,ex2].WQQ[k+p,0]


