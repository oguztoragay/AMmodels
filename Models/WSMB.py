# Updated on Feb 08 2021 for warm start for both linear and quadratic model, # Includes Nonlinear models (Baron) as well
# Includes all the figures writen in PDF format (all what we import from WSDW)
import time
import numpy as np
import pandas as pd
from openpyxl import Workbook
import GSgenerator as GS
from WSMILP import MILP, MILP_without
from Warm import warm
from WSQuadratic1 import GBNLPpyo, QUAD_without
from MINLP import NLPpyo
from WSDW import Draw_MINLP, Draw_MILP, Draw_Warm, Draw_GROUND_solid, Draw_GROUND_dashed
import os
if __name__ == '__main__':
    counter = 0
    try:
        r1_set = [[0, 0.2, 0.5], [0, 0.2, 0.3, 0.4, 0.5], [0, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]]
        r2_set = [0.2, 0.3, 0.5]
        E = 109000
        smax = 100000
        dmax = 0.095
        Load = [50]
        df = pd.DataFrame(columns=['Model', 'Load', 'LB', 'UB', 'Weight', 'Time', 'WS', 'Gap'])
        for ii in Load:
            print('<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')
            print('0000000000000000000000000000000000         WARM START ', str(ii), '           0000000000000000000000000000000000')
            print('<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')
            ins_f = (4, 4, [0, 3], [14], [ii], 1)
            # path_oguz = 'C:/Users/ozt0008/Documents/OneDrive - Auburn University/1 AM/Models/model PYTHON/8 February 2021/02.08.2021/results/'
            path_armin = 'C:/Users/ozt0008/OneDrive - Auburn University/1 AM/salam/'
            foldername = np.str(path_armin + str(ins_f[0]) + 'x' + str(ins_f[1]) + 'x' + str(ins_f[5]) + 'z')
#            foldername = 'C:/Users/ozt0008/Desktop/bashe'
            # foldername = 'C:/Users/ozt0008/OneDrive - Auburn University/1 AM/Paper 1/Lastplot'
            try:
                os.mkdir(foldername)
            except:
                foldername = str(foldername + '/' + time.strftime("%Y%m%d-%H%M%S"))
                os.mkdir(foldername)
            GS_ins_f = GS.Generate(ins_f[0], ins_f[1], ins_f[2], ins_f[3], ins_f[4], E, [0, 0.5])
            nodes = GS_ins_f.nodes
            elements = GS_ins_f.elements
            Xw, Yw, Dw, Vw, Sw, weight, TWS = warm(E, nodes, elements, [0, 0.5], dmax, smax)
            fname = np.str(str(ins_f[0]) + 'x' + str(ins_f[1]) + '_' + 'WSG' + '_' + str(ii))
            fname_ground = np.str(str(ins_f[0]) + str(ins_f[1]) + str(ins_f[5]))
            Draw_GROUND_dashed(nodes, elements, fname_ground, foldername)
            Draw_Warm(nodes, elements, Xw, TWS, weight, ii, fname, foldername)
            print('Warm Start took {} seconds.'.format(TWS))
            for jj in r1_set:
                print('<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')
                print('000000000000000000000000000000000         LINEAR WITH WS ', str(len(jj)), ' load: ', str(ii), '          00000000000000000000000000000000')
                print('<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')
                ins = (4, 4, [0, 3], [14], [ii], 1)
                GS_ins = GS.Generate(ins[0], ins[1], ins[2], ins[3], ins[4], E, jj)
                nodes = GS_ins.nodes
                elements = GS_ins.elements
                celements = GS_ins.celements
                localtime = time.asctime(time.localtime(time.time()))
                # ================================================================================
                print(':::: Current Time :::: {}'.format(localtime))
                print(':::: Current Instance is (with WS):{} :::: '.format(str(ins)))
                print(':::: Current Cross-sectional area set is {} :::: '.format(str(jj)))
                X, S, W, TLP, data1_wsmilp, data2_wsmilp = MILP(E, nodes, elements, jj, dmax, smax, Xw, Yw, Dw, Vw, Sw)
                fname = np.str(str(ins_f[0]) + 'x' + str(ins_f[1]) + '_' + 'WL' + '_' + str(len(jj)) + '_' + str(ii))
                Draw_MILP(nodes, elements, X, S, W, TLP, ii, jj, 'WS-MILP', dmax, fname, foldername)
                counter = counter + 1
                df.loc[counter] = ['CS' + str(r1_set.index(jj)+1), ii, data1_wsmilp[0].lower_bound, data1_wsmilp[0].upper_bound, W, data2_wsmilp[0].Time, 1, 0]
                # ================================================================================
                print('<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')
                print('0000000000000000000000000000000         LINEAR WITHOUT WS ', str(len(jj)), ' load: ', str(ii), '          0000000000000000000000000000000')
                print('<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')
                X1, S1, W1, TLP1, data1_milp, data2_milp = MILP_without(E, nodes, elements, jj, dmax, smax)
                fname = np.str(str(ins_f[0]) + 'x' + str(ins_f[1]) + '_' + 'L' + '_' + str(len(jj)) + '_' + str(ii))
                Draw_MILP(nodes, elements, X1, S1, W1, TLP1, ii, jj, 'MILP', dmax, fname, foldername)
                counter = counter + 1
                df.loc[counter] = ['CS' + str(r1_set.index(jj)+1), ii, data1_milp[0].lower_bound, data1_milp[0].upper_bound, W1, data2_milp[0].Time, 0, 0]
            # ============================================================================
            print('<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')
            print('000000000000000000000000000000          QUADRATIC MODEL_WS ', ' load: ', str(ii), '           000000000000000000000000000000')
            print('<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')
            Z, W, TQP, data_wsq = GBNLPpyo(E, nodes, celements, r2_set, dmax, smax, Xw, Yw, Dw, Vw, Sw)
            fname = np.str(str(ins_f[0]) + 'x' + str(ins_f[1]) + '_' + 'WQ' + '_' + str(ii))
            counter = counter + 1
            df.loc[counter] = ['WSQ', ii, data_wsq[0], data_wsq[1], data_wsq[2], data_wsq[4], 1, data_wsq[3]]
            Draw_MINLP(nodes, celements, Z, W, TQP, ii, 'QD', dmax, fname, foldername)
            # ================================================================================
            print('<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')
            print('0000000000000000000000000000000          QUADRATIC MODEL ', ' load: ', str(ii), '           0000000000000000000000000000000')
            print('<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')
            ZQ1, WQ1, TQP1, data_q = QUAD_without(E, nodes, celements, r2_set, dmax, smax)
            fname = np.str(str(ins_f[0]) + 'x' + str(ins_f[1]) + '_' + 'Q' + '_' + str(ii))
            counter = counter + 1
            df.loc[counter] = ['Q', ii, data_q[0], data_q[1], data_q[2], data_q[4], 0, data_q[3]]
            Draw_MINLP(nodes, celements, ZQ1, WQ1, TQP1, ii, 'QUAD', dmax, fname, foldername)
            # ================================================================================
            print('<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')
            print('0000000000000000000000000000000000       NONLINEAR MODEL ', ' load: ', str(ii), '        0000000000000000000000000000000000')
            print('<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')
            NLPstart = time.time()
            Z, W, TNLP, data1_minlp, data2_minlp = NLPpyo(E, nodes, celements, r2_set, dmax, smax, 'BARON', 'PC')
            fname = np.str(str(ins_f[0]) + 'x' + str(ins_f[1]) + '_' + 'NL' + '_' + str(ii))
            Draw_MINLP(nodes, celements, Z, W, TNLP, ii, 'NL', dmax, fname, foldername)
            counter = counter + 1
            df.loc[counter] = ['NL', ii, data1_minlp[0].lower_bound, data1_minlp[0].upper_bound, W, data2_minlp[0].Time, 0, 0]
        df.to_excel(foldername+'/'+str(ins_f[0]) + 'x' + str(ins_f[1]) + 'x' + str(ins_f[5])+'.xlsx')
    except ValueError as e:
        print(e)
# instances_3x3 = {1: (3, 3, [0, 2], [7], [ii], 1), 2: (3, 3, [0, 1, 3], [8], [ii], 2), 3: (3, 3, [0, 3, 6], [5], [-ii], 3)}
# instances_4x4 = {1: (4, 4, [0, 3], [14], [ii], 1), 2: (4, 4, [0, 1, 4], [15], [ii], 2), 3: (4, 4, [0, 4, 8, 12], [11], [-ii], 3)}
# instances_5x5 = {1: (5, 5, [0, 4], [22], [ii], 1), 2: (5, 5, [0, 1, 2, 5, 10], [24], [ii], 2), 3: (5, 5, [0, 5, 10, 15, 20], [14], [-ii], 3)}
