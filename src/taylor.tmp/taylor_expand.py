#!/usr/bin/python
# -*- python -*-   
# $RCSfile: taylor_expand.py,v $
"""Usage: python taylor_expand.py ncmExpand.txt

Expand and simplify input expressions to first/second order in the exteneral
momentum. Expressions are ported to FORM to handle the non-commutative products
of Dirac matrices easily.

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
def format_initial(content):
    #collect information from in_file
    str_expr = ""

    #format the input by removing the \'s at the end of each string
    #all the mathematica input stretches over multiple lines e.g.
    #  " expr1 \
    #    more of expr1 "
    #want to collect everything inbetween the ""'s and call that one string, 
    #and then append that string to one vector.
    for lines in range(0,len(content)):
        tmp_txt = content[lines]
        tmp_txt = tmp_txt.strip('\n') #remove any newlines
        tmp_txt = tmp_txt.strip(' \ ') #sometimes mathematica will place these to indicate EOL
	tmp_txt = tmp_txt.replace(" ", "") #remove all white space 
        tmp_txt = tmp_txt.replace(".","*")
        if tmp_txt.find("\"") != -1: 
            print "ERROR: wrong formatting in ncmExpand.txt, use ToString['expr',FormatType->InputForm]"
            exit()
	str_expr += tmp_txt #final expression is a string.	
    

    #collect all the lorentz indices in the numerator
    #this will include external indices as well
    #return the result as a dictionary, to be used later to append sums and contract repeated indices.

    lorentz_indices = collect_indices(str_expr,"[") #collects the lorentz indices

    #uncontract all repeated indices
    #str_expr = uncontract(str_expr, lorentz_indices)
    

    #call FORM to properly expand all the gamma matrices
    #WITHOUT simplifying them, yet
    #mathematica doesn't distribute non-commuting objects
    #easily, which is why FORM is used at this point

    str_expr = format_form(str_expr)
    str_expr = form_expand(str_expr, lorentz_indices)


    return str_expr

#----------------------------------------------------------------
def collect_indices(expr, delimit):
    #collects the lorentz indices from each monomial
    #and returns the result as a dictionary, 

    # arg : delimit is a str, assumes delimit = "[" for now
    # arg : monomial_vec = vector of strings
    
    lorentz_dict = {} #create the empty dictionary
    
    while expr.find(delimit) != -1: #finds all instances of "[" in monomial
        expr_new = expr.split(delimit) #splits the monomial at each instance of "["
        expr = expr.replace(delimit,"(")
        if len(expr_new) > 1: #if we have found some instances of "[ ]" in the monomial
            for args in expr_new:
                if args.find("]") != -1 and args.find(",") != -1: #if we can find an instance of "]" in the split string and , then we have del[ix1,ix2]
                    lorentz_key = args[0:args.find(",")]
                    lorentz_dict[lorentz_key] = 1
                    lorentz_key = args[args.find(",")+1:args.find("]")]
                    lorentz_dict[lorentz_key] = 1
                if args.find("]") != -1 and args.find(",") == -1:
                    lorentz_key = args[0:args.find("]")] # collects everything inside [ ... ] in monomial
                    
                    # remove any unwanted indices here.
                    if lorentz_key.find("+")!=-1 or lorentz_key.find("-")!=-1: 
                        continue
                    else:
                        lorentz_dict[lorentz_key] = 1 #append the result to a dictionary
    
    #transfer these to a vector for easier use
    #in the future and return all indices
    lorentz_indices = []
    for key in lorentz_dict:
    #    if key == 'y':
    #        print "ERROR: index y is reserved as an internal index in FORM"
    #        print "please choose a different index"
    #        exit()
        lorentz_indices.append(key)
    

    #quick fix to collect the lorentz indices properly
    #the F[a,b,c] structures cause problems..
    #lorentz_indices_new = []
    #for index in lorentz_indices:
    #    idx = index.split(',')
    #    for idx2 in idx:
    #        if len(idx2)>1:
    #            lorentz_indices_new.append(idx2)

    return lorentz_indices

#---------------------------------------------------------------
def format_form(str_arr):

    #make everything lowercase
    if len(str_arr) == 2:
        str = str_arr[1]
    else:
        str = str_arr
    str = str.replace("SIN", "sin")
    str = str.replace("Sin", "sin")
    str = str.replace("Cos", "cos")
    str = str.replace("COS", "cos")
    str = str.replace(" PLUS ", "+")
    str = str.replace(" MINUS ", "-")
    str = str.replace("MINUS","-")
    str = str.replace("I","i_")
    str = str.replace("i*","i_*")
    str = str.replace("(i)","(i_)")
    str = str.replace(".","*")
    str = str.replace("[","(")
    str = str.replace("]",")")
    str = str.replace(" ", "")
    
    if len(str_arr) == 2:
        str_formatted =[]
        str_formatted.append(str_arr[0])
        str_formatted.append(str)
    else:
        str_formatted = str

    return str_formatted

#----------------------------------------------------------------

def form_expand(input_str,lorentz_indices):

    # order of operations
    # 0. drop powers of 'a'
    # 1. sum over lorentz indices -- removing dirac deltas
    #    next step is to create a good sum routine for deltas.
    # 2. apply all possible trig identities -- removing cos for sin etc.
    # 3. power counting routine
    # 4. dirac algebra

    #construct the index list to be used:
    lorentz_str = ""
    for s in lorentz_indices:
        lorentz_str += s + "," #weird way of constructing a list of all lorentz indices
    if lorentz_str == "":
        print "ERROR: no lorentz indices found! Check mathematica.txt"
        exit()
    else:
        lorentz_str = lorentz_str[0:len(lorentz_str)-1] + ";\n"


    #simple code to make the file tmp.frm
    #used for FORM to make some final contractions/simplifications 
    input_str = input_str.replace(" ","")
    form_file = open('tmp.frm', 'w')
    form_file.write("off statistics;\n")
    form_file.write("autodeclare index dum,tmp;\n")
    form_file.write("Indices" + " " + lorentz_str)
    form_file.write("CFunctions p,k,kh,sin,cos,M,S,B,b,d,SUM,sum,I,F,f,sig,eps,s,m;\n")
    form_file.write("CFunctions J1,...,J1000;\n")
    form_file.write("Function G,g,T,t;\n")
    form_file.write("Symbol g0,AA,aa,rh,BB,bb,rw,m0,DB,db,DQ,dq,a,pow,sf;\n")
    form_file.write("Local str = "  + input_str + ";\n")
    form_file.write("repeat;\n")
    form_file.write("id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?)*G(dum5?)*G(dum6?)*G(dum7?)*G(dum8?) = G(dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8);\n")
    form_file.write("id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?)*G(dum5?)*G(dum6?)*G(dum7?) = G(dum1,dum2,dum3,dum4,dum5,dum6,dum7);\n")
    form_file.write("id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?)*G(dum5?)*G(dum6?) = G(dum1,dum2,dum3,dum4,dum5,dum6);\n")
    form_file.write("id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?)*G(dum5?) = G(dum1,dum2,dum3,dum4,dum5);\n")
    form_file.write("id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?) = G(dum1,dum2,dum3,dum4);\n")
    form_file.write("id G(dum1?)*G(dum2?)*G(dum3?) = G(dum1,dum2,dum3);\n")
    form_file.write("id G(dum1?)*G(dum2?) = G(dum1,dum2);\n")
    form_file.write("endrepeat;\n")

    form_file.write("contract;\n")
    form_file.write("print str;\n")
    form_file.write("*********\n")
    form_file.write("*ANSWER\n")
    form_file.write(".end\n")
    form_file.close()

    #now call form tmp.frm and save the output to be formatted.
    #NB: for some reason, mathematica does not like python to call FORM directly
    #    the I/O formatting between these two causes bugs, but developing a separate
    #    shell script tends to work.
    expr = os.popen("./form_expand.sh","r").read()

    #replace this comment if not using mathematica
    expr = expr.replace(" ", "")
    expr = expr.replace("\n","")
    expr = expr.replace("str=","")
    expr = expr.replace(";","")
    
    if expr == "":
        print "ERROR: form tmp.frm did not execute properly"
        expr = os.popen("form tmp.frm | grep Line","r").read()
        print "///"
        print expr
        exit()

    return expr

#----------------------------------------------------------------
def collect_args(str, fnc_str, delimit):

    # this routine actually sums over 
    # the delta functions, needs to be renamed.

    #str = expression
    # fncs = ['Sin','Cos',..etc.]
    # delimit = ["[","]"]
    #now separate each monomial in the numerator
    #for monomial in input_arr:
        #ugly string parsing ensues
        #first find all instances of sin and cos in monomial

    count =0
    place =0
    # I removed the while statement here since it was
    # slowing things way down for longer expressions
    # there are some new caveats with the speedier code:
    # sin/cos args are assumed to be 300 chars or less, 
    # although this is easily lengthened
    args=[]
    d_args = []

    for fncs in fnc_str:
        
        #find all places in str where fncs appear
        place = []
        for m in re.finditer(fncs, str):
            place.append( m.start())
        for strt in place:
            count = 0 
            substr = ""
 
            start= strt+len(fncs)
            tst = str[start:start+300]
            end = start + tst.find(")")+1
            
            arg = str[start:end]
            if arg.find(',')!=-1:
                arg = arg.strip('()')
                d_arg = arg.split(',')
                d_args.append(d_arg)
            else:
                args.append(arg.strip('()'))
            
                   
    # this is a much faster replacement routine 
    # than calling str.replace(x,y), which was VERY slow 
    # on large expressions...

    #str = mysubst(str,olds,reps) 

    return (args,d_args)

#---------------------------------------------------------------------------------

def mysubst(somestr, oldvec, subs):

    #create dictionary here
    vec = []
    for iter in range(0,len(oldvec)):
        vec.append((oldvec[iter],subs[iter]))

    dic = dict(vec)

    # use regex.
    pattern = re.compile('|'.join(re.escape(k) for k in dic))

    return pattern.sub(lambda k:dic[k.group()], somestr)

#---------------------------------------------------------------------

def vectorize(indata):

    # this routine is OK speed-wise

    formatted_data = []
    #separate terms beginning with +s and -s
    tmp_txt = indata.split('+')

    #formatted into tmp_txt = (-term1-term2, -term3, term4, term5-term6-term7-term8...etc.)
    tmp_txt_minus=[]
    for elmnt in range(len(tmp_txt)):
        if tmp_txt[elmnt].find('^-') != -1:
            tmp_txt[elmnt] = tmp_txt[elmnt].replace('^-', '{POW}')

        if tmp_txt[elmnt].find('-') != -1:
            #split e.g. term5-term6-term7 --> tmp_list = (term5, term6, term7)
            tmp_list = tmp_txt[elmnt].split('-')
            for str in range(1,len(tmp_list)): #skip the first element -- it is a plus sign always
                tmp_list[str] = tmp_list[str].replace('{POW}','^-')
                tmp_txt_minus.append(tmp_list[str])
                
                #replace in tmp_txt term5-term6-term7-.. ---> term5
                tmp_txt[elmnt] = tmp_list[0]


    tmp_txt_plus = []
    for elmnt in tmp_txt:
        if elmnt != "":
            elmnt = elmnt.replace('{POW}','^-')
            tmp_txt_plus.append(elmnt)
    #append the + signed data first
    #so formatted_data[0] = all terms starting with + sign
    formatted_data.append(tmp_txt_plus)
    #next append the - signed data
    #formatted_data[1] = all terms starting with - sign
    formatted_data.append(tmp_txt_minus)


    return formatted_data

#----------------------------------------------------------------
def dirac_collect(str):

#this routine will collect the dirac structure
#in each monomial and append the result to the 
#beginning of the string
   dirac_vec = []
   d_prod  = ""
   
   #routine to construct g_(1,x1,x2,x3..etc.) based on the list of G's appearing
   #in numerator expression
   while str.find("G(") != -1: #while there are some G's in this string
       str_new = str[str.find("G(")+2:] #collect the stuff inside G(...)
       str_new = str_new[:str_new.find(")")]
       dirac_vec.append(str_new)
       str = str.replace("G(" + str_new + ")","1") #WARNING: This may cause errors if G appears first in the char string


   for d_mtrx in range(0,len(dirac_vec)): #this code converts the G(x1).G(x2)... --> g_(1,x1,x2..)
       if d_mtrx != len(dirac_vec)-1:
           d_prod += dirac_vec[d_mtrx] + ","
       else:
           d_prod += dirac_vec[d_mtrx]
   d_prod = "g_(1," + d_prod + ");|"

   if d_prod == "g_(1,);|": d_prod = "g_(1,1);|" #simple error checking in case no dirac matrices present

   #make the final result
   dirac_str = d_prod + str

   return dirac_str

#----------------------------------------------------------------
def to_mathematica(input_arr):


    input_arr = input_arr.replace("(","[")
    input_arr = input_arr.replace(")","]")
    input_arr = input_arr.lower()
    input_arr = input_arr.replace("@","(")
    input_arr = input_arr.replace("??",")")
    input_arr = input_arr.replace("t[","T[")
    input_arr = input_arr.replace("f[","F[")
    input_arr = input_arr.replace("cd[","CD[")
    input_arr = input_arr.replace("cos","Cos")
    input_arr = input_arr.replace("sin","Sin")
    input_arr = input_arr.replace("plus","+")
    input_arr = input_arr.replace("minus","-")
    input_arr = input_arr.replace("i_","I")
    input_arr = input_arr.replace("g_[1,1]","1")
    input_arr = input_arr.replace("g_[1,","G[")
    input_arr = input_arr.replace("g_[1 ,","G[")
    input_arr = input_arr.replace("e_","eps")
    input_arr = input_arr.replace("d_","del")
    input_arr = input_arr.replace("5_","5")
    input_arr = input_arr.replace("denov","DenOv")


    #format the products of su(3) matrices
    #input_arr = color_format(input_arr)

    return input_arr
         
#----------------------------------------------------------------
#     BEGIN MAIN CODE
#----------------------------------------------------------------

if __name__ == "__main__":
 from sys import argv,exit

 
#-----pass options-------
# boolean var for summing over lorentz indices is the only option 
 parser = OptionParser()
 parser.add_option("-s", "--sum", action="store_true", dest="sum_opt", default = False,
                  help="sum over repeated Lorentz indices, leaving external indices untouched")

 parser.add_option("-f", "--file", action="store", dest="filename",default="ncmExpand.txt",
                  help="store a mathematica string to a text file to be read as input.")

 parser.add_option("-t", "--trace", action="store_true", dest="trace", default = False,
                   help="trace over lorentz indices in numerator.")

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
 
 
 expanded_form = format_initial(indata)
 mathematica_expanded = to_mathematica(expanded_form)


 print mathematica_expanded

# End $RCSfile: convert_form_new.py,v $
