# Updated on AUG 27
import time
import smtplib
import numpy as np
import GSgenerator as GS
from MILP import MILPpyo
#from for_gurobi import GBNLPpyo
#from MINLP import NLPpyo
from DW import Draw_GS, Draw_MILP, Draw_MINLP 
if __name__ == '__main__':
    try:
        r1_set = [0,0.3,0.4,0.5] # Set of possible radius for MILP model
        r2_set = [0.3, 0.5, 0.4] # Lower and Upper bound of radius in MINLP model [min radius, max redius, initial radius]
        E = 250000
        smax = 300
        dmax = 0.5 # Maximum displacement in the nodes connected to the structure
        Load = 2000
        ins = (5,3,[0,4],[12],[-Load])
        Problem_start_time = time.time()
        GS_ins = GS.Generate(ins[0],ins[1],ins[2],ins[3],ins[4],E,r1_set) # Structure to be solved
        nodes = GS_ins.nodes;     elements = GS_ins.elements;     celements = GS_ins.celements
        GS_finish_time = time.time() - Problem_start_time
        localtime = time.asctime( time.localtime(time.time()))
        print(':::: Current Time :::: {}'.format(localtime))
        print('Below Ground structure has been created. Total time for "GS" generating is "{}" seconds.'.format(np.round(GS_finish_time,2)))
        print(' 1 - Model to solve: {}\n 2 - Potential cross-sections set as: {}\n 3 - Nonlinear cs set: {}\n 4 - dmax = {} and smax = {} and load = {}'.format(ins, r1_set, r2_set, dmax, smax, Load))
        Draw_GS(nodes,elements)       
#======================================================================================================================================================================================
        print('----------------------------------------( LINEAR MODEL PC CPLEX )----------------------------------------')
        LPstart = time.time()
        #Options: cplex(NEOS), octeract-engine(PC), gurobi_ampl(PC), gurobi(PC)
        X, W = MILPpyo(E, nodes, elements, r1_set, dmax, smax, 'CPLEX', 'PC')# 
        TLP = np.round(time.time()-LPstart,3)
        Draw_MILP(nodes, elements, X, W, TLP, 'CP')
#======================================================================================================================================================================================
#        print('----------------------------------------( LINEAR MODEL PC GUROBI )----------------------------------------')
#        LPstart = time.time()
#        #Options: cplex(NEOS), octeract-engine(PC), gurobi_ampl(PC), gurobi(PC)
#        X, W = MILPpyo(E, nodes, elements, r1_set, dmax, smax, 'GUROBI', 'PC')# 
#        TLP = np.round(time.time()-LPstart,3)
#        Draw_MILP(nodes, elements, X, W, TLP, 'GB')
#======================================================================================================================================================================================
#        print('----------------------------------------( NONLINEAR MODEL PC GUROBI )-------------------------------------')
#        NLPstart = time.time()
#        #Options: BARON(PC), knitro(PC), knitro(NEOS), APOPT(APOPT), octeract-engine(PC)
#        Z, W = GBNLPpyo(E, nodes, celements, r2_set, dmax, smax, 'gurobi_ampl', 'PC')
#        TNLP = np.round(time.time()-NLPstart,3)
#        Draw_MINLP(nodes, celements, Z, W, TNLP, 'GU')
#======================================================================================================================================================================================            
#        print('----------------------------------------( NONLINEAR MODEL PC )--------------------------------------------')
#        NLPstart = time.time()
#        #Options: BARON(PC), knitro(PC), knitro(NEOS), APOPT(APOPT), octeract-engine(PC)
#        Z, W = NLPpyo(E, nodes, celements, r2_set, dmax, smax, 'BARON', 'PC')
#        TNLP = np.round(time.time()-NLPstart,3)
#        Draw_MINLP(nodes, celements, Z, W, TNLP, 'BA')
#======================================================================================================================================================================================
#        mail = smtplib.SMTP('smtp.gmail.com',587) 
#        mail.ehlo() 
#        mail.starttls() 
#        mail.login('oguztoragay@gmail.com','wkbipvwzrgdwdsbb') 
#        content = str('Model is done with all methods at: '+ str(time.asctime(time.localtime(time.time()))))
#        mail.sendmail('oguztoragay@gmail.com','ozt0008@auburn.edu', content)
#        mail.close()
#        print('*********************************************************************************************************')
#        print("Email was sent to ozt0008@auburn.edu to inform the user about the process of solving instance.")
#        print('*********************************************************************************************************')
    except ValueError as e:
        print(e)

# try LANCELOT on NEOS server. LANCELOT is suitable for large nonlinearly-constrained problems. 
# The basic algorithm combines the objective function and the set of all constraints more 
# complicated than simple bounds on the variables in an augmented Lagrangian function.
# A sequence of problems is solved, in which the current augmented Lagrangian function is
# approximately minimized within the region defined by the simple bounds.
  
 
