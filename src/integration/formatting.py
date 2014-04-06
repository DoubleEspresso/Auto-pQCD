#!/usr/bin/python

import fractions
#=================================================
# set of scripts to write/read from/to 
# FORM.  Used to port sub-expressions of a 
# lattice diagram to and from mathematica/VEGAS
# to be integrated.
#=================================================

class Parser:
    
    def __init__(self):  
        self.file = "tmp.out";

    def read(self,str):
        # simple script to format the 
        # form input, remove spaces etc.
        indata=[]
        expr=""
        if str == 'FORM':
            for line in open("tmp.out", 'rb'):        
                indata.append(line.strip().replace("\n",""))
                
            for line in indata:
                expr += line

        if str == 'MATHEMATICA':
            for line in open("math.in", 'rb'):
                indata.append(line.strip().replace("\n"),"")

            for line in indata:
                expr += line
                                
        return expr

    
    def vectorize(self,expr):
        # script to split the expression into 
        # signed monomial terms, both plus signed and 
        # minus signed data are saved in a vector and returned
        vector_data = []
        expr = expr.replace(' ','')
        expr = expr.replace('^-','{POW}')

        tmp = expr.split('+')

        tmp_minus = []
        tmp_plus  = []

        # negative terms
        for index in range(len(tmp)):
            if tmp[index].find('-') != -1:
                tmp_negative = tmp[index].split('-')
                for j in range(1,len(tmp_negative)):
                    tmp_negative[j] = tmp_negative[j].replace('{POW}','^-')
                    tmp_minus.append(tmp_negative[j])
                tmp[index] = tmp_negative[0] # the first element is always positive
                
        
        # positive terms
        for elmnt in tmp:
            if elmnt != "":
                elmnt = elmnt.replace('{POW}','^-')
                tmp_plus.append(elmnt)

        vector_data.append(tmp_plus)
        vector_data.append(tmp_minus)

        return vector_data


    def fetch_integrals(self,vector_data,options):
        # plus terms are vector_data[0]
        # mins terms are vector_data[1]
        B_INTS = []
        J_INTS = []
        F_INTS = []
        
        plus_vec = []
        mins_vec = []
        # this loop re-orders the powers appearing in F,J,B integrals
        # in ascending order and makes the replacements in the original expression
        # it then forms a list of unique integrations that must be done 
        # for overlap fermions we don't sort the aa,bb powers            
        for plus_monomial in vector_data[0]:
            if plus_monomial.find('B')!=-1:
                str = ""
                end_str = ""
                tmp1 = plus_monomial.split('B', 1)[1]
                tmp2 = tmp1.split(')',1)[0]
                str += tmp2[:4]
                if (len(tmp2)<=12): 
                    sorted_tmp2 = sorted(list(tmp2[4:].replace(",","")),key=int)
                if (len(tmp2)>12): 
                    sorted_tmp2 = sorted(list(tmp2[4:12].replace(",","")),key=int)
                    end_str = ","+tmp2[13:16]
                for pow in sorted_tmp2:
                    str += "," + pow
                B_INTS.append('B'+str + end_str+')')
                plus_vec.append(plus_monomial.replace('B'+tmp2+')','B'+str+end_str+')'))

            if plus_monomial.find('J')!=-1:
                str = ""    
                end_str = ""
                tmp1 = plus_monomial.split('J',1)[1]
                tmp2 = tmp1.split(')',1)[0]
                str += tmp2[:4]
                if (len(tmp2)<=12): 
                    sorted_tmp2 = sorted(list(tmp2[4:].replace(",","")),key=int)
                if (len(tmp2)>12): 
                    sorted_tmp2 = sorted(list(tmp2[4:12].replace(",","")),key=int)
                    end_str = ","+tmp2[13:16]
                for pow in sorted_tmp2:
                    str += "," + pow
                J_INTS.append('J'+str+end_str+')')
                plus_vec.append(plus_monomial.replace('J'+tmp2+')','J'+str+end_str+')'))
                
            if plus_monomial.find('F')!=-1:
                str=""
                end_str = ""
                tmp1 = plus_monomial.split('F',1)[1]
                tmp2 = tmp1.split(')',1)[0]
                str += tmp2[:4]
                if tmp2[4:].find("/") != -1:
                    expr=tmp2[4:].split(",")[1:]
                    sorted_tmp2=sorted(expr,key=eval)
                else:
                    if (len(tmp2)<=12): 
                        sorted_tmp2 = sorted(list(tmp2[4:].replace(",","")),key=int)
                        end_str = ""
                    if (len(tmp2)>12): 
                        sorted_tmp2 = sorted(list(tmp2[4:12].replace(",","")),key=int)
                        end_str = ","+tmp2[13:16]
                    for pow in sorted_tmp2:
                        str += "," + pow
                F_INTS.append('F'+str+end_str+')')
                plus_vec.append(plus_monomial.replace('F'+tmp2+')','F'+str+end_str+')'))

        for mins_monomial in vector_data[1]:
            if mins_monomial.find('B')!=-1:
                str = ""
                tmp1 = mins_monomial.split('B', 1)[1]
                tmp2 = tmp1.split(')',1)[0]
                str += tmp2[:4]
                if (len(tmp2)<=12): 
                    sorted_tmp2 = sorted(list(tmp2[4:].replace(",","")),key=int)
                    end_str = ""
                if (len(tmp2)>12): 
                    sorted_tmp2 = sorted(list(tmp2[4:12].replace(",","")),key=int)
                    end_str = ","+tmp2[13:16]
                for pow in sorted_tmp2:
                    str += "," + pow
                B_INTS.append('B'+str + end_str+')')
                mins_vec.append(mins_monomial.replace('B'+tmp2+')','B'+str+end_str+')'))

            if mins_monomial.find('J')!=-1:
                str = ""
                tmp1 = mins_monomial.split('J',1)[1]
                tmp2 = tmp1.split(')',1)[0]
                str += tmp2[:4]
                if tmp2[4:].find("/") != -1:
                    expr=tmp2[4:].split(",")[1:]
                    sorted_tmp2=sorted(expr,key=eval)
                else:
                    if (len(tmp2)<=12): 
                        sorted_tmp2 = sorted(list(tmp2[4:].replace(",","")),key=int)
                        end_str = ""
                    if (len(tmp2)>12): 
                        sorted_tmp2 = sorted(list(tmp2[4:12].replace(",","")),key=int)
                        end_str = ","+tmp2[13:16]
                    for pow in sorted_tmp2:
                        str += "," + pow
                    J_INTS.append('J'+str+end_str+')')
                    mins_vec.append(mins_monomial.replace('J'+tmp2+')','J'+str+end_str+')'))
                
            if mins_monomial.find('F')!=-1:
                str = ""
                tmp1 = mins_monomial.split('F',1)[1]
                tmp2 = tmp1.split(')',1)[0]
                str += tmp2[:4]
                if (len(tmp2)<=12): 
                    sorted_tmp2 = sorted(list(tmp2[4:].replace(",","")),key=int)
                    end_str = ""
                if (len(tmp2)>12): 
                    sorted_tmp2 = sorted(list(tmp2[4:12].replace(",","")),key=int)
                    end_str = ","+tmp2[13:16]
                for pow in sorted_tmp2:
                    str += "," + pow
                F_INTS.append('F'+str+end_str+')')
                mins_vec.append(mins_monomial.replace('F'+tmp2+')','F'+str+end_str+')'))

        # write the 'sorted' expression to file (mathematica formatting)
        self.write(plus_vec,mins_vec)
        #print '--B-integrals--'
        #print list(set(B_INTS))
        #print ''
        #print '--J-integrals--'
        #print list(set(J_INTS))
        #print ''
        #print '--F-integrals--'
        #print list(set(F_INTS))
        #exit()
        # remove duplicates using set.
        return (list(set(B_INTS)), list(set(J_INTS)), list(set(F_INTS)))
        


    def write(self, plus_terms, minus_terms):
        # format all monomial terms to mathematica 
        # recombine the plus/minus data into one string
        # and write the result to file

        final_str = ""
        for monomial in plus_terms:
            monomial = monomial.replace("(","[")
            monomial = monomial.replace(")","]")
            monomial = monomial.replace("i_","I")
            final_str += monomial + "+"

        for monomial in minus_terms:
            monomial = monomial.replace("(","[")
            monomial = monomial.replace(")","]")
            monomial = monomial.replace("i_","I")
            final_str += "-" + monomial
        
        final_str = final_str.replace("+-","-")
        

        OUT_FILE = open("out/lattice_integrand.out","w")
        OUT_FILE.write(final_str)
        OUT_FILE.close()
        return
            
                
    
    def export(self,integral_vector):
        B_INTS = integral_vector[0]
        J_INTS = integral_vector[1]
        F_INTS = integral_vector[2]

        B_FILE = open("out/b_integrals.out", "w")
        J_FILE = open("out/j_integrals.out", "w")
        F_FILE = open("out/f_integrals.out", "w")

        for int in B_INTS:
            B_FILE.write(int+'\n')
        for int in J_INTS:
            J_FILE.write(int+'\n')
        for int in F_INTS:
            F_FILE.write(int+'\n')

        B_FILE.close()
        J_FILE.close()
        F_FILE.close()

        return




    def convert(self,expr,type):
        # short script to format mathematica code
        # to all FORM conventions.
        expr = expr.replace("[","(")
        expr = expr.replace("]",")")
        return
        
