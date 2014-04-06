#!/usr/bin/python
# -*- python -*-   
# $RCSfile: main.py,v $
"""
Usage: python main.py --options

  Options:
    --trace     : will compute traces over all gamma matrices present in the expression.
    --integrate : will computing all integrals, default behavior is no numerical integration.    

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

#-----------------------------------------
import integration/formatting



#------------------------------------------
#- Main -
#------------------------------------------
    
if __name__ == "__main__":
 from sys import argv,exit

#-----pass options-------
 parser = OptionParser()
 
 parser.add_option("-t", "--trace", action="store", dest="filename",default="Formatted_Data.txt",
                   help="stores resulting formatted data in given file.")
                   
                   
 (options, args) = parser.parse_args()


#---------------------------------------------------
#-- parse the header file, set external indices etc.


#-- format the mathematica input for form.
Formatter = formatting.Parser()
mathData  = Formatter.read('MATHEMATICA')
formData  = Formatter.convert(mathData,'FORM')
Formatter.write_to_form()

# run form for various options.
