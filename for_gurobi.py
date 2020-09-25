# 09/08/2020
import warnings
warnings.filterwarnings("ignore", category = UserWarning)
import numpy as np
import gurobipy as gb
from gurobipy import GRB
from shapely.geometry import LineString
import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
# %% Definig Sets, Variables and Parameters
def GBNLPpyo(E, nodes, celements,r2_set, dmax, smax, sol, wheresol):
    boundary = []; load_node = []; nfree = []; load_node_dof = []; load_value = []
    for i in nodes.keys():
        if nodes[i].tip == 1:
            boundary.append(i)
            nfree.append(nodes[i].dof)
        elif nodes[i].tip == 2:
            load_node.append(i)
            load_node_dof.append(nodes[i].dof[1])
            load_value.append(nodes[i].load)
    f = np.ravel([nodes[i].f for i in nodes.keys()])
## --------------------------------------------------------------------------------------------------------------------------------------------------
    m       = gb.Model('NewModel')
    m.setParam('NonConvex', 2)
    LN    = list(nodes.keys());    LE    = list(celements.keys())
    dofs  = list(np.ravel([nodes[i].dof for i in nodes.keys()]));    nfree = list(np.ravel(nfree));    free  = list(filter(lambda i: i not in nfree, dofs))
    amin  = np.pi*(r2_set[0]**2);  amax  = np.pi*(r2_set[1]**2)
## --------------------------------------------------------------------------------------------------------------------------------------------------
    d   = m.addVars(dofs, lb = -dmax, ub = dmax, vtype = GRB.CONTINUOUS, name = 'Disp')
    for i in nfree:
        d[i].ub = 0
        d[i].lb = 0
    ar  = m.addVars(LE, lb = amin, ub = amax, vtype = GRB.SEMICONT, name = 'Area')
    ar_2= m.addVars(LE, vtype = GRB.CONTINUOUS, name = 'Area^2')
    RF = m.addVars(nfree, vtype = GRB.CONTINUOUS)
    le  = {}
    for i in LE:
        le.update({i:celements[i].length}) 
# %% Generating Objective function and Constraints
    obj = gb.quicksum(ar[i]*le[i] for i in LE)
    m.setObjective(obj, GRB.MINIMIZE)
#-------------------------------------------------------------CONSTRAINTS---------------------------------------------------------------------------      
    for satr in free: # equilibrium equations
        temp1 = 0; temp2 = 0; temp3 = 0
        for i in LE:
            temp1 += celements[i].KE[0]*ar[i]  *celements[i].B[0][satr]*gb.quicksum(celements[i].B[0][di0]*d[di0] for di0 in free)
            temp2 += celements[i].KE[1]*ar_2[i]*celements[i].B[1][satr]*gb.quicksum(celements[i].B[1][di1]*d[di1] for di1 in free)
            temp3 += celements[i].KE[2]*ar_2[i]*celements[i].B[2][satr]*gb.quicksum(celements[i].B[2][di2]*d[di2] for di2 in free)
        m.addConstr(temp1 + temp2 + temp3 == f[satr], name='Equilibrium')
#-----------------------------------------------------------------------------------------
    for i in LE:
        m.addQConstr(ar_2[i] == ar[i]*ar[i], name='PowerOfAreas')
#-----------------------------------------------------------------------------------------      
    temp1 = [x for x in LN if x not in boundary]
    LNf = [x for x in temp1 if x not in load_node]
    for i in LNf: # displacement on each DoF limited by the existance of at least one beam
        dofs = nodes[i].dof
        ww = nodes[i].where
        for ii in dofs:
            m.addConstr((d[ii]<=sum(ar[j] for j in ww)*5), name='dispULimit')
            m.addConstr((d[ii]>=-sum(ar[j] for j in ww)*5), name='disLLimit')
#-----------------------------------------------------------------------------------------   
    for i in LE: # Removing the cross-overing beams
        for j in range(i+1, len(LE)):
            seg1 = LineString([[celements[i].nodei.x, celements[i].nodei.y], [celements[i].nodej.x, celements[i].nodej.y]])
            seg2 = LineString([[celements[j].nodei.x, celements[j].nodei.y], [celements[j].nodej.x, celements[j].nodej.y]])
            int_pt = seg1.intersects(seg2)
            toucher = seg1.touches(seg2)
            if int_pt == True and toucher == False:
               m.addConstr(ar[i]*ar[j] == 0, name='Crossing')
#-----------------------------------------------------------------------------------------               
#    for i in LE: # Stress bounds tried...
#        dofi = nodes[celements[i].orient[0]].dof
#        dofj = nodes[celements[i].orient[1]].dof
#        c = celements[i].cosan
#        s = celements[i].sinan
#        di = c*d[dofi[0]]+s*d[dofi[1]]
#        dj = c*d[dofj[0]]+s*d[dofj[1]]   
#        m.addConstr((dj-di)*E<=smax*ar[i])#(celements[i].KE[0])
#        m.addConstr(-smax*ar[i]<=(dj-di)*E)           
    m.update() 
#    m.display()
# %% Solving MINLP model    
#    m.linearize()
#    m.setParam(GRB.Param.Threads, 20.0)
    m.setParam(GRB.Param.TimeLimit, 36000)
    m.optimize()
    weight1 = m.objVal
    force_loc = {}
    for i in LE:
        dofz = celements[i].dof
        d1 = d[dofz[0]].x;        d2 = d[dofz[1]].x
        d3 = d[dofz[2]].x;        d4 = d[dofz[3]].x
        d5 = d[dofz[4]].x;        d6 = d[dofz[5]].x
        dz = [d1,d2,d3,d4,d5,d6]
        KEL = np.multiply(celements[i].ke1*(ar[i].x),celements[i].bloc[0])+ \
              np.multiply(celements[i].ke2*(ar[i].x)**2,celements[i].bloc[1])+ \
              np.multiply(celements[i].ke3*(ar[i].x)**2,celements[i].bloc[2])
        KELT = np.dot(celements[i].transform.T,np.dot(KEL,celements[i].transform))
        F = np.dot(KELT,dz)
        force_loc[i] = F.T
    print('-------------------------------------------( BEAMS )-------------------------------------------')
    Y = []
    for i in LE:
        Y.append(ar[i].x)
    toprint = pd.DataFrame(index = LE, columns = ['No','Area','Radius', 'N1','V1','M1','N2','V2','M2'])
    for i in LE:
        toprint.loc[i] = [i,ar[i].x,np.sqrt(ar[i].x/np.pi),force_loc[i][0],force_loc[i][1],force_loc[i][2],force_loc[i][3],force_loc[i][4],force_loc[i][5]] 
    indexNames = toprint[(toprint['Area'] == 0)].index
    toprint.drop(indexNames, inplace=True)
    print(toprint)
    print('-------------------------------------------( NODES )-------------------------------------------')
    toprint2 = pd.DataFrame(index = nodes.keys(), columns = ['No','U1','U2','UR3'])
    for i in nodes.keys():
        dofz = nodes[i].dof
        dis = [d[dofz[0]].x, d[dofz[1]].x, d[dofz[2]].x]
        toprint2.loc[i] = [i,dis[0],dis[1],dis[2]]
    print(toprint2)
    print('-------------------------------------------( REACTIONS )---------------------------------------')    
    RF = {}
    for satr in nfree:
        temp11 = 0; temp21 = 0; temp31 = 0
        for i in LE:
            temp11 += celements[i].KE[0]*ar[i].x  *celements[i].B[0][satr]*gb.quicksum(celements[i].B[0][di0]*d[di0].x for di0 in dofs)
            temp21 += celements[i].KE[1]*ar_2[i].x*celements[i].B[1][satr]*gb.quicksum(celements[i].B[1][di1]*d[di1].x for di1 in dofs)
            temp31 += celements[i].KE[2]*ar_2[i].x*celements[i].B[2][satr]*gb.quicksum(celements[i].B[2][di2]*d[di2].x for di2 in dofs)
        RF[satr]= temp11 + temp21 + temp31
    for nn in boundary:
        doffs = nodes[nn].dof
        print('RF in node {} is: in X direction: "{}", in Y direction: "{}" and rotation: "{}"'.format(nn,RF[doffs[0]],RF[doffs[1]],RF[doffs[2]]))   
#    print(Y)
    return(Y,weight1)