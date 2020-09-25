# Updated on SEP 12 for DoE
import time
import numpy as np
import GSgenerator as GS
from MILP import MILPpyo
from for_gurobi import GBNLPpyo
from MINLP import NLPpyo
from DWdetail import Draw_MILP, Draw_MINLP 
if __name__ == '__main__':
    try:
        r1_set = [[0,0.2,0.5], [0,0.2,0.3,0.4,0.5], [0,0.2,0.25,0.3,0.35,0.4,0.45,0.5]]
        r2_set = [0.2, 0.5, 0.4]
        E = 250000
        smax = 100000
        dmax = 0.25
        Load = [1000]#,250,500]
        for ii in Load:
            for jj in r1_set:
                ins = (3,3,[0,2],[7],[ii])
#                ins = (3,3,[0,1,3],[8],[ii])
#                ins = (3,3,[0,3,6],[5],[-ii])
                GS_ins = GS.Generate(ins[0],ins[1],ins[2],ins[3],ins[4],E,jj)
                nodes = GS_ins.nodes     
                elements = GS_ins.elements
                celements = GS_ins.celements
                localtime = time.asctime( time.localtime(time.time()))
                print(':::: Current Time :::: {}'.format(localtime))
                LPstart = time.time()
                X, W = MILPpyo(E, nodes, elements, jj, dmax, smax, 'GUROBI', 'PC')
                TLP = np.round(time.time()-LPstart,3)
                Draw_MILP(nodes, elements, X, W, TLP, ii, jj,'GB',dmax)
#            #================================================================================
            QPstart = time.time()
            Z, W = GBNLPpyo(E, nodes, celements, r2_set, dmax, smax, 'GUROBI', 'PC')
            TNLP = np.round(time.time()-QPstart,3)
            Draw_MINLP(nodes, celements, Z, W, TNLP, ii, 'QG', dmax)
            #================================================================================
#            NLPstart = time.time()
#            Z, W = NLPpyo(E, nodes, celements, r2_set, dmax, smax, 'BARON', 'PC')
#            TNLP = np.round(time.time()-NLPstart,3)
#            Draw_MINLP(nodes, celements, Z, W, TNLP, ii, 'BA', dmax)
    except ValueError as e:
        print(e)