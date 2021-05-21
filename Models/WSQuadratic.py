# 12/07/2020 Updated to use the warm start for quadratic model
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import numpy as np
import gurobipy as gb
from gurobipy import GRB
from shapely.geometry import LineString
import pandas as pd
import time

def GBNLPpyo(E, nodes, celements, r2_set, dmax, smax, Xw, Yw, Dw, Vw, Sw):
    boundary = []
    load_node = []
    nfree = []
    load_node_dof = []
    load_value = []
    for i in nodes.keys():
        if nodes[i].tip == 1:
            boundary.append(i)
            nfree.append(nodes[i].dof)
        elif nodes[i].tip == 2:
            load_node.append(i)
            load_node_dof.append(nodes[i].dof[1])
            load_value.append(nodes[i].load)
    f = np.ravel([nodes[i].f for i in nodes.keys()])
    # --------------------------------------------------------------------------------------------------------------------------------------------------
    m = gb.Model('NewModel')
    m.setParam('NonConvex', 2)
    m.setParam('TimeLimit', 18000)  # Time limit increase from 10800 to 18000 on 11.23.2020
    m.setParam('Threads', 16)  # Input your machine's # of CPUs
    m.setParam('LogToConsole', 1)
    m.setParam('DisplayInterval', 100)
    m.setParam('InfUnbdInfo', 1)
    m.setParam('FeasibilityTol', 1e-5)
    m.setParam('Heuristics', 1)
    m.setParam('RINS', 10)
    m.setParam('Cuts', 3)
    m.setParam('MIPFocus', 1)
    m.setParam('Presolve', 2)
    m.setParam('PreQLinearize', 2)
    # --------------------------------------------------------------------------------------------------------------------------------------------------
    LN = list(nodes.keys())
    LE = list(celements.keys())
    dofs = list(np.ravel([nodes[i].dof for i in nodes.keys()]))
    nfree = list(np.ravel(nfree))
    free = list(filter(lambda i: i not in nfree, dofs))
    amin = np.pi * (r2_set[0] ** 2)
    amax = np.pi * (r2_set[2] ** 2)
    # --------------------------------------------------------------------------------------------------------------------------------------------------
    d = m.addVars(dofs, lb=-dmax, ub=dmax, vtype=GRB.CONTINUOUS, name='Disp')
    for i in nfree:
        d[i].ub = 0
        d[i].lb = 0
    ar = m.addVars(LE, lb=amin, ub=amax, vtype=GRB.SEMICONT, name='Area')
    ar_2 = m.addVars(LE, vtype=GRB.CONTINUOUS, name='Area^2')
    y = m.addVars(LN, vtype=GRB.BINARY, name='Node_appear')
    le = {}
    for i in LE:
        le.update({i: celements[i].length})
    obj = gb.quicksum(ar[i] * le[i] for i in LE)
    m.setObjective(obj, GRB.MINIMIZE)
    # -------------------------------------------------------------CONSTRAINTS---------------------------------------------------------------------------
    for satr in free:  # equilibrium equations
        ps = nodes[np.floor(satr / 3)].where
        temp1 = 0
        temp2 = 0
        temp3 = 0
        for i in ps:  # LE
            temp1 += celements[i].KE[0] * ar[i] * celements[i].B[0][satr] * gb.quicksum(celements[i].B[0][di0] * d[di0] for di0 in free)
            temp2 += celements[i].KE[1] * ar_2[i] * celements[i].B[1][satr] * gb.quicksum(celements[i].B[1][di1] * d[di1] for di1 in free)
            temp3 += celements[i].KE[2] * ar_2[i] * celements[i].B[2][satr] * gb.quicksum(celements[i].B[2][di2] * d[di2] for di2 in free)
        m.addQConstr(temp1 + temp2 + temp3 <= f[satr], name='Equilibrium1'+str(satr))
        m.addQConstr(temp1 + temp2 + temp3 >= f[satr], name='Equilibrium2'+str(satr))
    # -----------------------------------------------------------------------------------------
    for i in LE:
        m.addQConstr(ar_2[i] <= ar[i] * ar[i], name='PowerOfAreas1'+str(i))
        m.addQConstr(ar_2[i] >= ar[i] * ar[i], name='PowerOfAreas2'+str(i))
    # -----------------------------------------------------------------------------------------
    LNf = set(LN) - set(boundary) - set(load_node)
    for i in LNf:  # displacement on each DoF limited by the existence of at least one beam
        dofs = nodes[i].dof
        ww = nodes[i].where
        for ii in dofs:
            m.addConstr((d[ii] <= sum(ar[j] for j in ww)), name='dispULimit'+str(ii))
            m.addConstr((d[ii] >= -sum(ar[j] for j in ww)), name='dispLLimit'+str(ii))
    # -----------------------------------------------------------------------------------------
    for i in LE:  # Removing the cross-overing beams
        for j in range(i + 1, len(LE)):
            seg1 = LineString(
                [[celements[i].nodei.x, celements[i].nodei.y], [celements[i].nodej.x, celements[i].nodej.y]])
            seg2 = LineString(
                [[celements[j].nodei.x, celements[j].nodei.y], [celements[j].nodej.x, celements[j].nodej.y]])
            int_pt = seg1.intersects(seg2)
            toucher = seg1.touches(seg2)
            if int_pt == True and toucher == False:
                m.addQConstr(ar[i] * ar[j] <= 0, name='Crossing1['+str(i)+','+str(j)+']')
                m.addQConstr(ar[i] * ar[j] >= 0, name='Crossing2['+str(i)+','+str(j)+']')
    # -----------------------------------------------------------------------------------------            
    for i in LE:
        m.addQConstr(ar[i] - amax * y[celements[i].nodei.name] * y[celements[i].nodej.name], rhs=0, sense=GRB.LESS_EQUAL, name='salam1'+str(i))
    # ----------------------------------------------------------
    for i in LNf:
        elemlist = nodes[i].where
        m.addConstr(sum(ar[elem] for elem in elemlist) >= 2 * y[i] * amin, name='salam2'+str(i))
    # ----------------------------------------------------------
    m.update()
    for i in LE:
        if Xw[i, 1] == 1:
            ar[i].setAttr('Start', amax)
            ar_2[i].setAttr('Start', amax**2)
        else:
            ar[i].setAttr('Start', 0)
            ar_2[i].setAttr('Start', 0)
    for j in LN:
        if Yw[j] == 1:
            y[j].setAttr('Start', 1)
        else:
            y[j].setAttr('Start', 0)
    m.update()
    QPstart = time.time()
    m.optimize()
    TQP = np.round(time.time() - QPstart, 3)
    m.write("result1.mst")
    data_wsq = [m.ObjBound, m.ObjBoundC, m.ObjVal, m.MIPGap, m.Runtime]

    weight1 = m.objVal
    force_loc = {}
    for i in LE:
        dofz = celements[i].dof
        d1 = d[dofz[0]].x
        d2 = d[dofz[1]].x
        d3 = d[dofz[2]].x
        d4 = d[dofz[3]].x
        d5 = d[dofz[4]].x
        d6 = d[dofz[5]].x
        dz = [d1, d2, d3, d4, d5, d6]
        KEL = np.multiply(celements[i].ke1 * ar[i].x, celements[i].bloc[0]) + np.multiply(
            celements[i].ke2 * ar[i].x ** 2, celements[i].bloc[1]) + np.multiply(celements[i].ke3 * ar[i].x ** 2,
                                                                                 celements[i].bloc[2])
        KELT = np.dot(celements[i].transform.T, np.dot(KEL, celements[i].transform))
        F = np.dot(KELT, dz)
        force_loc[i] = F.T

    navadanam = np.zeros(len(LE))
    stress = np.zeros(len(LE))
    for i in LE:
        strain = 0
        if ar[i].x != 0:
            dof1 = celements[i].dof
            dal = []
            for j in dof1:
                dal.append(d[j].x)
            lenprime = np.sqrt(((dal[3] + celements[i].nodej.x) - (dal[0] + celements[i].nodei.x)) ** 2 + (
                    (dal[4] + celements[i].nodej.y) - (dal[1] + celements[i].nodei.y)) ** 2)
            strain = (lenprime - celements[i].length) / celements[i].length
            stress[i] = strain * E
            navadanam[i] = stress[i] * ar[i].x
    print('-------------------------------------------( BEAMS )-------------------------------------------')
    Y = []
    nodey = []
    for i in LE:
        Y.append(ar[i].x)
    for j in LN:
        nodey.append(y[j].x)
    toprint = pd.DataFrame(index=LE, columns=['No', 'Area', 'Radius', 'N1', 'V1', 'M1', 'N2', 'V2', 'M2', 'navadanam', 'Stress'])
    for i in LE:
        toprint.loc[i] = [i, ar[i].x, np.sqrt(ar[i].x / np.pi), force_loc[i][0], force_loc[i][1], force_loc[i][2],
                          force_loc[i][3], force_loc[i][4], force_loc[i][5], navadanam[i], stress[i]]
    indexNames = toprint[(toprint['Area'] == 0)].index
    toprint.drop(indexNames, inplace=True)
    print(toprint)
    print('-------------------------------------------( NODES )-------------------------------------------')
    toprint2 = pd.DataFrame(index=nodes.keys(), columns=['No', 'U1', 'U2', 'UR3'])
    for i in nodes.keys():
        dofz = nodes[i].dof
        dis = [d[dofz[0]].x, d[dofz[1]].x, d[dofz[2]].x]
        toprint2.loc[i] = [i, dis[0], dis[1], dis[2]]
    print(toprint2)
    print(Y)
    print('-----------------------------------------------------------------------------------------------')
    return Y, weight1, TQP, data_wsq