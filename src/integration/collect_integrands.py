#!/usr/bin/python
# -*- python -*-   
# $RCSfile: collect_integrands.py,v $
"""
Usage: python collect_integrands integrals.dat

Create some replacement rules for quick input to mathematica.
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

#-----------------------------------
import formatting
import integrate

#------------------------------------------
#- Main -
#------------------------------------------
    
if __name__ == "__collect_integrands__":
 from sys import argv,exit

#-----pass options-------
 parser = OptionParser()
 
 parser.add_option("-f", "--file", action="store", dest="filename",default="Formatted_Data.txt",
                   help="stores resulting formatted data in given file.")
                   
                   
 (options, args) = parser.parse_args()

#--------------------------------------------------
#-- format the input from FORM
Formatter = formatting.Parser()
formData = Formatter.read('FORM') 
vectorData = Formatter.vectorize(formData)

#-- get the integrals
integralVector = Formatter.fetch_integrals(vectorData)

#-- print them to files b_integrals.out, j_integrals.out
#-- include options(Action) here for overlap case!
Formatter.export(integralVector)



#-- chose a method to numerically compute the j-integrals
Integrator = integrate.Integrator()
j_integrals = Integrator.mathematica_integrate()


# -- use form to compute the b-integrals analytically
b_integrals = Integrator.form_integrate()

Integrator.export(j_integrals+b_integrals)

# -- prep the integrals.results file for mathematica replacements--
Integrator.format_results()

# -- run mathematica to replace all integrals, expand in epsilon, and 
#    print the result to final.result.out
ans = Integrator.mathematica_expand()


# print the final answer
print ""
print "=========================="
print "Done!"
print "=========================="
print ""
print "J = " + ans
print ""
exit()


