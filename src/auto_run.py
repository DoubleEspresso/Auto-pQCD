#!/usr/bin/python
# -*- python -*-   
# $RCSfile: main.py,v $
"""
Usage: python auto_run.py [Options]

Wrapper script to main.py for a range of Wilson rw and rho
parameters.  We can scan over rw in 0.1 increments for instance
and write all results to a single file in one go, rather than 
changing the rw parameter manually for each iteration.

Options:
  --f,file  : writes the result to user defined file.
  
"""

_author__ = "M. Glatzmaier"
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 1/1/2012 16:00:22 $"
__copyright__ = "Copyright (C) 2012 M. Glatzmaier"
__licence__ = """This file is covered by the GPL v2 (see LICENSE file)."""

import sys
import os
import re
import math
import subprocess
from subprocess import Popen
from optparse import OptionParser

if __name__=="__main__":
    #from sys import argv, exit
#----pass options----
    parser = OptionParser()
    parser.add_option("--file", action="store", dest="RESULT_FILE",default="auto_run_results.txt",help="stores all results of the auto run.")
    parser.add_option("--log", action="store", dest="LOG_FILE",default="auto_run.log",help="log file for the python code.")
    parser.add_option("--rw",type="float",action="store",dest="START_RW",help="set the initial value of the Wilson rw parameter.")
    parser.add_option("--max_rw",type="float",action="store",dest="MAX_RW",help="set the max value of the Wilson rw parameter.")
    parser.add_option("--delta_rw",type="float",action="store",dest="DELTA_RW",help="set the change in rw for each iteration.")
    parser.add_option("--rho",type="float",action="store",dest="START_RHO",help="set the initial value of the Overlap rho parameter.")
    parser.add_option("--max_rho",type="float",action="store",dest="MAX_RHO",help="set the max value of the Overlap rho parameter.")
    parser.add_option("--delta_rho",type="float",action="store",dest="DELTA_RHO",help="set the change in rho for each iteration.")
    parser.add_option("--only_integrate",type="int",action="store",dest="ONLY_INTEGRATE",help="specify if we are only integrating previously calculated results")

# -- on user input error -- 
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()    
    
    (options, args) = parser.parse_args()
    

# relabel the options args
START_RW       = options.START_RW
MAX_RW         = options.MAX_RW
DELTA_RW       = options.DELTA_RW
RESULT_FILE    = options.RESULT_FILE
LOG_FILE       = options.LOG_FILE
ONLY_INTEGRATE = options.ONLY_INTEGRATE

if (RESULT_FILE == "") :
    print "ERROR :: must pass a file to write results!"
    exit()

if (ONLY_INTEGRATE == ""):
    print "ERROR :: must specify if integrating only."
    exit()

# define the number of iterations
ITERATIONS = int((MAX_RW - START_RW)/DELTA_RW)+1


DIAGRAM_IN_HEADER= """# define the header diagram.in file
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
$$ Action        ::  Wilson
$$ TaylorOrder   ::  1
$$ Method        ::  mathematica
$$ Integrate     ::  true
$$ Error         ::  6
$$ Trace         ::  true
"""

DIAGRAM_IN_FOOTER="""
$$ rho           ::  0.5
$$ Sum           ::  
$$ ExtIndices    ::  
$$ ExtMomentum   ::  p
$$ Factor        ::  1/a^4
$$ Diagram       ::
#$$ Diagram       ::  WQQG[2*k+p,ex1].WQQ[k,0].QOP1[k].WQQ[k,0].WQQG[2*k+p,ex2].WQQ[k+p,0]"""



for iter in range(0,ITERATIONS):
    thisRW = START_RW + iter*DELTA_RW
    
    # we have the correct rw, we write this to 
    # diagram.in file and execute the code
    DIAGRAM_IN_STRING = "$$ rw            ::  " + str(thisRW)
    DIAGRAM_IN = DIAGRAM_IN_HEADER + DIAGRAM_IN_STRING + DIAGRAM_IN_FOOTER

    # write to diagram.in
    OUT_FILE = open("diagram.in","w")
    OUT_FILE.write(DIAGRAM_IN)
    OUT_FILE.close()
    
    # run the main file, dump all output to auto_run.log
    MAIN_EXE = ""
    if (ONLY_INTEGRATE == 0):
        MAIN_EXE = "./main.py"
        # only compute the integrand once,
        # use the results for each other iteration
        # ONLY_INTEGRATE = 1

    #if (ONLY_INTEGRATE == 1):
    #    MAIN_EXE = "./auto_integrate.py"
    if (MAIN_EXE == "") :
        print " ERROR :: incorrect setting for integrate options."
        exit()
    subprocess.call([MAIN_EXE])

    

    # we have returned from the calculation and now can store the 
    # results 
    ITER_DATA = " lattice integrand, iter :: "+ str(iter) + ", rw = " + str(thisRW) + " ,rho = --" + " twist_2_q->g sails diagram"
    # grab the result from the log file, and remove the log file
    with open(LOG_FILE) as RF:
        answer_lines = RF.readlines()
    RF.close()

    ANSWER = ITER_DATA+"\n"+answer_lines[len(answer_lines)-1]+"\n"+"\n"

    # append this result to the results file
    with open(RESULT_FILE,"a") as RESULTS:
        RESULTS.write(ANSWER)
    RESULTS.close()

#-------- END MAIN----------------    
    

