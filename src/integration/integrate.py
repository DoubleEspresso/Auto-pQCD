#!/usr/bin/python

import os
import subprocess
import threading
import multiprocessing
import time
#=====================================================
# set of integration scripts to read from 
# j_integrals.out, port them to mathematica/VEGAS
# and numerically integrate.  Results are formatted 
# for input into mathematica as replacement rules in 
# the final expression
#======================================================

# directories -- path to form executable
FORM_EXEC = "/home/mjg/physics/utils/form"

#----------------------------------------------------------------------
class myThread (threading.Thread):
    def __init__(self, threadID, j_ints, method, opts):
        threading.Thread.__init__(self)
        self.ID  = threadID
        self.j_ints = j_ints
        self.options=opts
        self.method =method

    def run(self):

        if self.method == 'mathematica':
            for j_int in self.j_ints:
                Integrator().integrate(j_int, "mathematica", self.ID, self.options)

        elif self.method == 'cubature':
            for j_int in self.j_ints:
                Integrator().integrate(j_int,"cubature",self.ID, self.options)

#----------------------------------------------------------------------

class Integrator:
    def __init__(self):
        self.file = "j_results.out"
        
    def mathematica_integrate(self,options):
        # we read in the J-integrals
        # from j_integrals.out
        indata = []
        j_answers = []
        iter = 1
        for line in open("out/j_integrals.out",'rb'):
            indata.append(line.strip().replace("\n",""))
            
        for line in open("out/f_integrals.out",'rb'):
            indata.append(line.strip().replace("\n",""))

        # execute parallel code here 
        nCPUS = multiprocessing.cpu_count()
        nINTS = len(indata)
        print ''
        print ' ...Integrating', nINTS, 'integrals'
        print ' ...using', nCPUS, 'CPUS.'
        print ''
        delta = int(nINTS / nCPUS)
        remainder = nINTS - nCPUS * delta
        lims = []
        threads = []
        if (nCPUS > 1):
            for cpu in range(1, nCPUS):
                lims.append( cpu * delta)
            lims.append(lims[len(lims)-1] + delta + remainder)
        else: threads.append(myThread(1,indata,options["Method"],options))

        if (nCPUS > 1):
            for thread in range(0,nCPUS):
                if thread == 0 : low = 0
                else : low = lims[thread-1]
                threads.append(myThread(thread,indata[low:lims[thread]],options["Method"],options))

        start_time = time.time()
        for thread in threads:
            thread.start()

        for thread in threads:
            thread.join()
        
        #------collect-------
        J_ANSWERS = []
        with open("math.txt") as f:
            content=f.readlines()
        for integral in content:
            J_ANSWERS.append(integral.replace("\"","").replace(" ","").replace("\n",""))
            
        subprocess.call(["rm","math.txt"])

        print ''
        print "Done."
        print "Process took: ", (time.time() - start_time)/60.0, "min."
        print ''

        return J_ANSWERS


    def form_integrate(self,options):
        # here we read the B-integrals
        # from b_integrals.out
        indata = []
        b_answers = []
        iter = 1

        for line in open("out/b_integrals.out",'rb'):
            indata.append(line.strip().replace("\n",""))

        # include coefficients for overlap integration here
        print indata
        print ''
        
        for b_integral in indata:
            print ''
            print 'Reducing:', b_integral
            print 'Integral:', iter, '/', len(indata)
            # format the integral for FORM, we order the arguments
            # such that the B(0,i,n1,n2,n3,n4) where n1<=n2<=n3<=n4
            # this operation is invariant under the hypercubic group and makes
            # the reduction algorithm in FORM easier to implement.
            first_piece = b_integral[:5]
            last_piece  = ""
            index_list = b_integral[6:].replace(")","")
            index_list = sorted(list(index_list.replace(",","")))
            for str in index_list:
                last_piece += str + ","
            last_piece = last_piece[:len(last_piece)-1] + ")"
            form_b = first_piece+","+last_piece
            iter+=1
            # only main thread executes these integrals, they are very
            # fast in general.
            b_answers.append(self.integrate(form_b,'form',1,options))

        return b_answers


    def integrate(self,j_int,method,ID,options):
        # the main numerical method, this is an internal method
        # used to call mathematica or cubature/VEGAS to compute the integral
        # numerically.

        # collect the args
        if (j_int.split("(")[0] == 'F') : fj_arg = "0"
        else:  fj_arg = "1"

        tmp = j_int.split("(")[1].strip(")")
        int_args = tmp.split(",")
        rw_arg   = options["rw"]
        rho_arg  = options["rho"]
        action   = options["Action"]
        id_arg   = str(ID)
        err_arg  = options["Error"]#"7"
        dq_arg   = int_args[0]
        dg_arg   = int_args[1]
        nx_arg   = int_args[2]
        ny_arg   = int_args[3]
        nz_arg   = int_args[4]
        nt_arg   = int_args[5]

        ## MATHEMATICA ROUTINE ##        
        if method == 'mathematica':

            if action  == 'Wilson' : 
                MATH_EXE = "./wilson_integrand.m"            
                subprocess.call([MATH_EXE,err_arg,dq_arg,dg_arg,nx_arg,ny_arg,nz_arg,nt_arg,rw_arg,rho_arg,fj_arg,id_arg])
            if action == 'Overlap' :
                MATH_EXE = "./overlap_integrand.m"
                if (len(int_args) >5): apow_arg = int_args[6]
                if (len(int_args) >6): bpow_arg = int_args[7]
                # if overlap action and J integral, append a coefficient arg to the integral
                # this coefficient is then computed in mathematica such that the numerical integration
                # is finite in the k->0 limit.  NB J is the difference in 2 integrals, so we set
                # the coefficient to ensure this difference is finite as k->0.
                #if fj_arg=="1":                            
                #    print "append coefficient to this integral!"
                subprocess.call([MATH_EXE,err_arg,dq_arg,dg_arg,nx_arg,ny_arg,nz_arg,nt_arg,rw_arg,rho_arg,fj_arg,id_arg,apow_arg,bpow_arg])
            # execute

            return

        elif method == "cubature":
            if (j_int.split("(")[0] == 'F') : f_or_j = 0.
            else:  f_or_j = 1.
            tmp = j_int.split("(")[1].strip(")")

            j_args = tmp.split(",")
            f_args = ''
            for j_arg in j_args:
                f_args += j_arg + " "

            tol = '1e-4'
            dim = '4'

            arg_str = dim + ' ' + tol + ' ' + f_args  + options["rw"] + " " + str(f_or_j)

            J_ANSWERS.append(j_int.replace("(","[").replace(")","]") +"->"+
                             os.popen("./wilson_integrand " + arg_str, "r").read().replace("\n","").replace("\11.","").replace("\"","")+",")
            return

        ## FORM INTEGRATION METHOD ##
        elif method == 'form':
            outfile = open('boson_reduction.frm')
            file_lines = outfile.readlines()
            outfile.close()
            if len(file_lines) == 0:
                print "Error : did not parse boson_reduction.frm correctly"
                print ""
                exit()
            
            for line in range(0,len(file_lines)):
                file_lines[line] = file_lines[line].strip('\n')
                
                # write the boson_integral to file
                if file_lines[line] == "***THE INTEGRAND FROM PYTHON***":
                    file_lines[line+1] = "" # erase this line
                    str1 = "L expr = " +j_int+ ";\n"
                    file_lines[line+1] = str1
                    str1 = ""

            # write to file
            outfile = open ('boson_reduction.frm','w')
            for line in range(0,len(file_lines)):
                print >> outfile, file_lines[line]
            outfile.close()

            # call form to compute b_integral analytically
            # this prints the result to bos.tmp, may want a 
            # shell script instead, so form ouput is surpressed.
            #call([FORM_EXEC, "boson_reduction.frm"])
            ans = os.popen("./form_run.sh","r").read()

            # format this answer for mathematica
            j_int = j_int.replace("(","[")
            j_int = j_int.replace(")","]")
            return j_int + "->" + ans.replace("\n","") + ","

        else:
            print "NOT IMPLEMENTED"
            exit()



    def export(self,results):
        INTS = open("out/integrals.results","w")

        for answer in results:
            INTS.write(answer+'\n')

        return


    def format_results(self):
        integrals = open("out/integrals.results")
        file_lines = integrals.readlines()
        integrals.close()

        for line in range(len(file_lines)):
            if line == 0: 
                str1 = "intReps = { " + file_lines[line]
                file_lines[line] = str1

            if line == len(file_lines)-1 : 
                str2 = file_lines[line].replace("\n","") + "};\n"
                file_lines[line] = str2.replace(",}","}")

        

        outfile = open("out/integrals.results","w")
        for line in range(0,len(file_lines)):
            print >> outfile, file_lines[line].strip("\n")
        outfile.close()
        return


    def mathematica_expand(self,options):
        # the final routine that replaces all integrals just computed
        # and expands everything in epsilon
        math_file = open("expand_integrand.m")
        file_lines = math_file.readlines()
        math_file.close()
        
        int_file = open("out/integrals.results")
        int_lines = int_file.readlines()
        int_file.close()

        integrals = ""
        for line in int_lines:
            integrals += line.strip("\n")
        
        num_file  = open("out/lattice_integrand.out")
        num_lines = num_file.readlines()
        num_file.close()

        integrand = ""
        for line in num_lines:
            integrand += line.strip("\n")

        for line in range(len(file_lines)):
            if line == 1:
                file_lines[line] = integrals
            
            if line == 2:
                file_lines[line] = "latticeIntegrand = " + integrand + ";\n"

        outfile = open("expand_integrand.m","w")
        for line in range(0,len(file_lines)):
            print >> outfile, file_lines[line].strip("\n")
        outfile.close()

        # set rw, rho 
        rw_arg   = options["rw"]
        rho_arg  = options["rho"]

        # run the file.
        MATH_EXE = "./expand_integrand.m"
        subprocess.call([MATH_EXE,rw_arg,rho_arg])
        with open("math.txt") as f:
            content = f.readlines()
        
        ans=""
        for stuff in content: 
            ans += stuff.replace(" ","").replace("\"","")

        subprocess.call(["rm","math.txt"])

        return ans
    #ans = os.popen("./math_run.sh expand","r").read()
    #    ans = ans.replace("\"","")

     #   return ans
