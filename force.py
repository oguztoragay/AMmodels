import time
import numpy as np
import GSgenerator as GS
from MILP import MILPpyo
from for_gurobi import GBNLPpyo
from MINLP import NLPpyo
from DW import Draw_GS, Draw_MILP, Draw_MINLP 
if __name__ == '__main__':
    try:
        r1_set = [0,0.3,0.4,0.5] # Set of possible radius for MILP model
        r2_set = [0.3, 0.5, 0.4] # Lower and Upper bound of radius in MINLP model [min radius, max redius, initial radius]
        E = 250000
        smax = 300
        dmax = 0.5 # Maximum displacement in the nodes connected to the structure
#        Load = 100
#        ins_dictionary = {0:(3,2,[0,2],[4],[Load]) ,1:(3,3,[0,2],[7],[Load]), 2:(3,3,[0,3,6],[2],[-Load]), 3:(3,3,[0,3,6],[7,8],[Load,Load]), 4:(3,3,[0,3,6],[5],[-Load]), 
#                          5:(3,4,[0,3,6,9],[5,8],[-Load,-Load]), 6:(3,4,[0,3,6,9],[5,10],[-Load,-Load]), 7:(4,3,[0,4,8],[7],[-Load]), 8:(4,3,[0,4,8],[3],[-Load]),
#                          9:(4,4,[0,4,8,12],[15],[Load]), 10:(4,4,[3,7],[13,14],[Load,Load]), 11:(4,4,[12,13],[3,7],[-Load,-Load]), 12:(4,4,[8,12],[2,3],[-Load,-Load]),
#                          13:(4,4,[12,13],[2,3],[-Load,-Load]), 14:(5,5,[0,1],[24],[Load])}
#        ins_dictionary ={2:(3,3,[0,3,6],[2],[-Load]), 3:(3,3,[0,3,6],[7,8],[Load,Load]), 7:(4,3,[0,4,8],[7],[-Load]), 9:(4,4,[0,4,8,12],[15],[Load])}
        for i in range(9):
#            ins_index = i        
#            ins = ins_dictionary[ins_index]
            Load = 100 + i*100
            ins = (3,3,[0,3,6],[5],[-Load])
            Problem_start_time = time.time()
            GS_ins = GS.Generate(ins[0],ins[1],ins[2],ins[3],ins[4],E,r1_set) # Structure to be solved
            nodes = GS_ins.nodes;     elements = GS_ins.elements;     celements = GS_ins.celements
            GS_finish_time = time.time() - Problem_start_time
            localtime = time.asctime( time.localtime(time.time()))
            print(':::: Current Time :::: {}'.format(localtime))
            print('Below Ground structure has been created for models number {}. Total time for "GS" generating is "{}" seconds.'.format(i,np.round(GS_finish_time,2)))
            print(' 1 - Model to solve: {}\n 2 - Potential cross-sections set as: {}\n 3 - Nonlinear cs set: {}\n 4 - dmax = {} and smax = {}'.format(ins, r1_set, r2_set, dmax, smax))
            Draw_GS(nodes,elements)       
#======================================================================================================================================================================================
#            print('----------------------------------------( LINEAR MODEL CPLEX )----------------------------------------')
#            LPstart = time.time()
#            #Options: cplex(NEOS), octeract-engine(PC), gurobi_ampl(PC), gurobi(PC)
#            X, W = MILPpyo(E, nodes, elements, r1_set, dmax, smax, 'CPLEX', 'PC')# 
#            TLP = np.round(time.time()-LPstart,3)
#            Draw_MILP(nodes, elements, X, W, TLP, 'CP')
#======================================================================================================================================================================================
            print('----------------------------------------( LINEAR MODEL GUROBI )----------------------------------------')
            LPstart = time.time()
            #Options: cplex(NEOS), octeract-engine(PC), gurobi_ampl(PC), gurobi(PC)
            X, W = MILPpyo(E, nodes, elements, r1_set, dmax, smax, 'gurobi', 'PC')# 
            TLP = np.round(time.time()-LPstart,3)
            Draw_MILP(nodes, elements, X, W, TLP, 'GB')
#======================================================================================================================================================================================
#            print('----------------------------------------( LINEAR MODEL XPRESS )----------------------------------------')
#            LPstart = time.time()
#            #Options: cplex(NEOS), octeract-engine(PC), gurobi_ampl(PC), gurobi(PC)
#            X, W = MILPpyo(E, nodes, elements, r1_set, dmax, smax, 'XPRESS', 'PC')# 
#            TLP = np.round(time.time()-LPstart,3)
#            Draw_MILP(nodes, elements, X, W, TLP, 'XP')
#======================================================================================================================================================================================
            print('----------------------------------------( NONLINEAR MODEL GUROBI )-------------------------------------')
#            NLPstart = time.time()
#            #Options: BARON(PC), knitro(PC), knitro(NEOS), APOPT(APOPT), octeract-engine(PC)
#            Z, W = GBNLPpyo(E, nodes, celements, r2_set, dmax, smax, 'gurobi_ampl', 'PC')
#            TNLP = np.round(time.time()-NLPstart,3)
#            Draw_MINLP(nodes, celements, Z, W, TNLP, 'GB')
#======================================================================================================================================================================================            
#            print('----------------------------------------( NONLINEAR MODEL BARON )--------------------------------------------')
#            NLPstart = time.time()
#            #Options: BARON(PC), knitro(PC), knitro(NEOS), APOPT(APOPT), octeract-engine(PC)
#            Z, W = NLPpyo(E, nodes, celements, r2_set, dmax, smax, 'BARON', 'PC')
#            TNLP = np.round(time.time()-NLPstart,3)
#            Draw_MINLP(nodes, celements, Z, W, TNLP, 'BA')
#======================================================================================================================================================================================
            print('************************************************************************************************************************')
    except ValueError as e:
        print(e)