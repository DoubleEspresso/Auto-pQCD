#!/usr/bin/python
# -*- python -*-   
# $RCSfile: io.py,v $

#---------------------------
# this is a library for 
# index routines:
#  1. index contraction
#  2. summing over kh-indices
#  3. misc tensor methods
#---------------------------
import re
import io

#-------------------------------------------
class Tensor:
    
    def __init__(self):
        self.file = "tmp.out"

#-------------------------------------------
    def contract(self, expr):

        IO = io.IO()

        # assumes expression is formatted for FORM 
        tmp = re.sub(r'ix(\d+)', r'%ix\1%',expr)
        vector_tmp = IO.vectorize(tmp)


        # the signed-monomials
        plus_contract = []
        mins_contract = []
        for iter in range(0,len(vector_tmp)):
            for monomial in vector_tmp[iter]:
                sum_inds   = self.f7(self.collect_indices(monomial,"SUM","FORM"))

                for sum_ind in sum_inds:
                    delta_inds = self.collect_indices(monomial,"DELTA","FORM")


                    for del_pair in delta_inds:
                        if (sum_ind == del_pair[0] or sum_ind == del_pair[1]):

                            # several cases
                            if (del_pair[0] == del_pair[1]) : 
                                if (monomial.rfind(del_pair[0]) == 3 and monomial.find('sum(')!=-1):
                                    monomial=monomial.replace('sum('+sum_ind+')','1')
                                    monomial=monomial.replace('d('+del_pair[0]+','+del_pair[0]+')','DIM')
                                else:
                                    monomial=monomial.replace("d("+del_pair[0]+","+del_pair[1]+")","1")
                            else:
                                if (sum_ind == del_pair[0]) : 
                                    monomial=monomial.replace("sum("+sum_ind+")","1",1)
                                    monomial=monomial.replace(del_pair[0],del_pair[1])
                                    monomial=monomial.replace("d("+ del_pair[1]+","+del_pair[1]+")","1")


                                    
                                if (sum_ind == del_pair[1]) : 
                                    monomial=monomial.replace("sum("+sum_ind+")","1",1)
                                    monomial=monomial.replace(sum_ind,del_pair[0])
                                    monomial=monomial.replace("d("+del_pair[0]+","+del_pair[0]+")","1")


                if (iter == 0) : plus_contract.append(monomial.replace("%",""))
                if (iter == 1) : mins_contract.append(monomial.replace("%",""))
     
        # check if the expression is non-zero
        if (len(plus_contract)==0 and len(mins_contract) == 0):
            print " Integrand vanished after Dirac simplification."
            exit()

        expr = IO.recombine(plus_contract,mins_contract)

        return expr

#-------------------------------------------
    def collect_indices(self, expr, type, format):
        
        indices = []
        
        if (format=="FORM"):        delimiters = ["(",")"]
        if (format=="MATHEMATICA"): delimiters = ["[","]"]

        # get the function type here
        if (type=="SUM"):
            fnc = 'sum'
        if (type=="DELTA"):
            fnc = 'd'
        if (type=="KH"):
            fnc = 'kh'
        if (type=="TRIG"):
            fnc = ['s','m','b','G']
        if (type=="DIRAC"):
            fnc = 'G'
        if (type=="EXTERNAL"):
            expr = re.sub(r'ex(\d+)',r'$ex\1$',expr)
            tmp = expr.split('$')[1:]
            for split_fncs in tmp:
                indices.append(split_fncs.split(delimiters[1])[0])
            return indices
            

        # append the indices by splitting the expression based on fnc
        if (type != "TRIG"):
            tmp = expr.split( fnc + delimiters[0])[1:]
            for split_fncs in tmp:
                if (type == "DELTA") : indices.append(split_fncs.split(delimiters[1])[0].split(','))
                else : indices.append(split_fncs.split(delimiters[1])[0])
                
        else:
            for f in fnc:
                tmp = expr.split( f + delimiters[0])[1:]
                for split_fncs in tmp:
                    indices.append(split_fncs.split(delimiters[1])[0])
            

        return indices
#-------------------------------------------
    def append_sums(self,expr,options):

        # want to vectorize the expression
        # loop over monomials
        # count occurence of sum_index
        # split off indices based on occurence
        # append sums
        #
        # pass to form to expand, recollect those sums, IR-reduce -- faster to do this one term at a time
        # for super large expressions form will hang
        #
        # recombine each monomial
        # return expression

        IO = io.IO()

        # assumes expression is formatted for FORM 
        vector_tmp = IO.vectorize(expr)

        # the signed-monomials
        plus_sums = []
        mins_sums = []
        for iter in range(0,len(vector_tmp)):

            if (iter == 0) : print "...Expanding sums for positive monomials"
            if (iter == 1) : print "...Expanding sums for negative monomials"

            for monomial in vector_tmp[iter]:                
                #print 'before replace:', monomial

                # this assumes expr in FORM format
                # get the kh-indices and sum-indices
                kh_inds = self.f7(self.collect_indices(monomial,"KH","FORM"))
                sum_inds= self.f7(self.collect_indices(monomial,"SUM","FORM"))


                # collect the powers of kh-here
                monomial = re.sub(r'\^(\d+)',r'^?\1?',monomial)
                monomial = re.sub(r'\^-(\d+)',r'^?-\1?',monomial)
        
                # append powers to a pow vector
                tmp  = monomial.split('^')[1:]
                pows = []

                for split_fncs in tmp:
                    if split_fncs.find('?')!=-1:
                        pows.append(split_fncs.split('?')[1])


                # re-label the kh's based on their powers
                # append sums here, for those monomials
                # that have a kh-index appearing in multiple places,
                # we remove the sum-over kh, giving kh a tmp index.
                for sum_index in sum_inds:
                    
                    for kh_index in kh_inds:
                        if (sum_index == kh_index):
                            # check if sum_index appears more than once in kh-fnc 
                            # if so, we give it a tmp index not equal to the external
                            # indices in the problem.
                            if ( monomial.count(sum_index) > 2 and monomial.find("kh("+sum_index+")")!= -1) :
                                for pow in pows:
                                    monomial = monomial.replace("kh("+sum_index+")^?"+pow+"?","kh(4)^?"+pow+"?",1)
                                break
                                
                            monomial = monomial.replace("sum("+sum_index+")","1",1)
                            for pow in pows:
                                monomial = monomial.replace("kh("+sum_index+")^?"+pow+"?","sum_("+sum_index+",1,4,kh( "+sum_index+ ")^"+pow+")",1)
                
                            # includes powers of 1
                            monomial = monomial.replace("kh("+sum_index+")","sum_("+sum_index+",1,4,kh( "+sum_index+ ")^"+pow+ ")",1)
                            if (monomial.find('kh( '+sum_index+"))^")!=-1): 
                                print 'ERROR: The term:'
                                print monomial
                                print 'is incorrect'
                                exit()

                # note on indices below : two classes of diagrams,
                # one including a G operator and one class of Q operators
                # for q->q and g->q there are no 'polarization indices' since there are no external gauge fields
                # but for q->g and g->g we have additional polarization indices of the gauge fields.  For these 
                # sets of diagrams we re-interpret the indices below to be external gauge field indices (polarization indices)
                # and set ex1=ex2=1 for forward matrix elements.  This amounts to searching for a del.p^2 tensor structure in the final
                # calculations for diagrams mixing into the G operator.

                # do the same for the external indices
                # need better code for this one....
                ex_inds = ['ex1','ex2','ex3','ex4']
                for ex_index in ex_inds:
                    if (ex_index == 'ex1') : monomial=monomial.replace(ex_index,'3')
                    if (ex_index == 'ex2') : monomial=monomial.replace(ex_index,'3')
                    if (ex_index == 'ex3') : monomial=monomial.replace(ex_index,'3')
                    if (ex_index == 'ex4') : monomial=monomial.replace(ex_index,'4')

                monomial = monomial.replace("?","")
                #if monomial != '0': print 'after replace:', monomial
                #if monomial != '0': print ''                

                # call form on this monomial
                if (options["Action"]=="Wilson"):
                    monomial = IO.write_form_file(monomial,"IR_ISOLATE")                
                if (options["Action"]=="Overlap"):
                    monomial = IO.write_form_file(monomial,"IR_ISOLATE_OVERLAP")
                if (options["Action"]!="Wilson" and options["Action"]!="Overlap"):
                    print " Error : unknown action in diagram.in, cannot proceed."
                    exit()
                if (monomial == '') : monomial = '0'
                #if monomial != '0': print 'after form:', monomial
                #if monomial != '0': print ''                


                if (iter == 0) : plus_sums.append(monomial)
                if (iter == 1) : mins_sums.append(monomial)

        expr = IO.recombine(plus_sums,mins_sums)


        return expr

#-------------------------------------------
# returns unique elements in a list
    def f7(self,seq):
        seen = set()
        seen_add = seen.add
        return [ x for x in seq if x not in seen and not seen_add(x)]
