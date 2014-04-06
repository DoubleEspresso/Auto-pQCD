#!/usr/bin/python
# -*- python -*-   
# $RCSfile: convert_form_new.py,v $
"""Usage: python subscript.py

Short script to append subscripts to aribtrary expressions
in a correct manner.  Mathematica does a surprisingly poor job
with this.

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

#----------------------------------------------------------------
#				DEFINITIONS
#----------------------------------------------------------------
def subscript(expr,script):
    # a character string of elements that should not
    # be subscripted.
    no_sub='123456789/*+-!()'

    res = ''
    for ch in expr:
        sub = 1
        for ns in no_sub:
            if ch == ns:
               sub = 0
        if sub == 1:
            ch = 'Subscript[' + ch + ',' + script + ']'
        res += ch
    
    return res

#----------------------------------------------------------------
#     BEGIN MAIN CODE
#----------------------------------------------------------------
if __name__ == "__main__":
 from sys import argv,exit

# possible options for file input
 parser = OptionParser()

 parser.add_option("-f", "--file", action="store", dest="filename",default="index.txt",
                  help="store a mathematica string to a text file to be read as input.")

 parser.add_option("-s","--sub",action="store",dest="sub_sc",default="NOINDEX",
                   help="The string to be used as a subscript in expression.")

 (options, args) = parser.parse_args()

# Parse options from command line:
# No options given -> print usage summary and exit
 if len(argv) > 4:
  print __doc__
  exit()
 
 # read input from mathematica
 in_file = open(options.filename)
 indata = in_file.readlines()
 in_file.close()


 # pass the mathematica input to subscript fnc
 # to correctly subscript everything.
 indata = indata[0].replace(" ","")
 print subscript(indata.strip('\n'),options.sub_sc)
 exit()
