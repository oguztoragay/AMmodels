# 06/24/2020 Most updated MILP
import warnings
warnings.filterwarnings("ignore", category = UserWarning)
import numpy as np
from shapely.geometry import LineString
from pyomo.environ import Param, ConcreteModel, Var, Objective, ConstraintList, value, minimize, Binary, Constraint
from pyomo.opt import SolverFactory
from pyomo.opt.parallel import SolverManagerFactory
# import winsound
from pyomo.util.infeasible import log_infeasible_constraints
import pandas as pd
import itertools
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
# %% Definig Sets, Variables and Parameters
def MILPpyo(E, nodes, elements, r_set, dmax, smax, sol, wheresol):
    boundary = set(); load_node = set(); nfree = set(); load_node_dof = []; load_value = []
    for i in nodes.keys():
        if nodes[i].tip == 1:
            boundary.add(i)
            nfree.update(nodes[i].dof)
        elif nodes[i].tip == 2:
            load_node.add(i)
            load_node_dof.append(nodes[i].dof[1])
            load_value.append(nodes[i].load)
    f = np.ravel([nodes[i].f for i in nodes.keys()])
#----------------------------------------------------------
    m       = ConcreteModel()
    m.LN    = set(nodes.keys()) # List of nodes
    m.LE    = set(elements.keys()) # List of elements
    m.PR    = set(range(len(r_set))) # List of profiles
    m.dofs  = set(np.ravel([nodes[i].dof for i in nodes.keys()]))
    m.nfree = nfree
    m.free  = m.dofs - m.nfree
    m.rho = 1; m.E = E; m.dmax = dmax; m.Smax = smax
#----------------------------------------------------------
    m.d  = Var(m.dofs, initialize = 0) # displacement along each degrees of freedom
    for i in m.nfree:
        m.d[i].fix(0)
    m.x   = Var(m.LE, m.PR, domain = Binary, initialize = 0) # index i is used for the elements or beams
    m.y   = Var(m.LN, domain = Binary, initialize = 0) # indicator of each nodes to be connected to the structure
    m.v   = Var(m.LE, {0,1,2}, m.PR, initialize = 0) # Elongation or contraction of beam i with profile p
    m.s   = Var(m.LE, {0,1,2}, initialize = 0) # forces (axial, shear or rotational) on each bar
    m.Q   = Param(m.LE, initialize = 0, mutable = True) # forces (axial or rotational) on each degree of freedom
    m.RF  = Var(m.nfree, initialize = 0) # Reaction forces on the boundary condition fixed degrees of freedom
    m.M   = Param(initialize = 1000) # Big M for the purpose of linearization of logic constraint
# --------------------------------------------------------------OBJECTIVE---------------------------------------------------------------------------
    m.obj = Objective(expr = sum([elements[i].profile[p].vol*m.x[i,p] for i in m.LE for p in m.PR]), sense = minimize) # Objective function minimizing the volum
#-------------------------------------------------------------CONSTRAINTS---------------------------------------------------------------------------    
    m.cons1 = ConstraintList()
    for satr in m.dofs:
        temp1 = 0
        for i in m.LE:
            for p in m.PR:
                temp1 += (elements[i].profile[p].KE[0])*(m.v[i,0,p])*(elements[i].B[0][satr]) + \
                         (elements[i].profile[p].KE[1])*(m.v[i,1,p])*(elements[i].B[1][satr]) + \
                         (elements[i].profile[p].KE[2])*(m.v[i,2,p])*(elements[i].B[2][satr])
        if satr in m.nfree:
            m.cons1.add(temp1 == m.RF[satr])
#            m.cons1.add(temp1 >= m.RF[satr])
#            Constraint.Skip
        else:
            m.cons1.add(temp1 == f[satr])
#            m.cons1.add(temp1 >= f[satr])
#----------------------------------------------------------
    m.cons2 = ConstraintList()
    for i in m.LE:
        for j in {0,1,2}:
            m.cons2.add(sum(m.v[i,j,p] for p in m.PR-{0})-sum(elements[i].B[j][d]*m.d[d] for d in m.dofs) <= m.M*(1-sum(m.x[i,p] for p in m.PR)))
            m.cons2.add(sum(m.v[i,j,p] for p in m.PR-{0})-sum(elements[i].B[j][d]*m.d[d] for d in m.dofs) >= -1*m.M*(1-sum(m.x[i,p] for p in m.PR)))
#----------------------------------------------------------    
    m.cons3 = ConstraintList()
    for i in m.LE:
        for j in {0,1,2}:
            m.cons3.add(elements[i].profile[p].KE[j]*(sum(m.v[i,j,p] for p in m.PR-{0})) - m.s[i,j] <= m.M*(1-sum(m.x[i,p] for p in m.PR-{0})))
            m.cons3.add(elements[i].profile[p].KE[j]*(sum(m.v[i,j,p] for p in m.PR-{0})) - m.s[i,j] >= -1*m.M*(1-sum(m.x[i,p] for p in m.PR-{0})))
##----------------------------------------------------------             
    m.cons4 = ConstraintList()
    for i in m.LE:
        for j in {0,1,2}:
            m.cons4.add(m.s[i,j]<= m.M*(sum(elements[i].profile[p].area*m.x[i,p] for p in m.PR-{0})))
            m.cons4.add(m.s[i,j]>= -m.M*(sum(elements[i].profile[p].area*m.x[i,p] for p in m.PR-{0})))
#----------------------------------------------------------             
    m.cons5 = ConstraintList()
    for i in m.LE:
        for p in m.PR-{0}:
            for j in {0,1,2}:
                m.cons5.add(m.v[i,j,p] <= m.M*m.x[i,p])
                m.cons5.add(m.v[i,j,p] >= -1*m.M*m.x[i,p])
#----------------------------------------------------------
    m.cons6 = ConstraintList()
    for i in m.LE:
        m.cons6.add(sum(m.x[i,p] for p in m.PR) <= 1)
#----------------------------------------------------------
    m.cons7 = ConstraintList()
    for i in m.LE:
        m.cons7.add(sum(m.x[i,p] for p in m.PR-{0}) <= 0.5*(m.y[elements[i].nodei.name] + m.y[elements[i].nodej.name]))
#----------------------------------------------------------        
    m.cons8 = ConstraintList()
    for i in m.LN-boundary-load_node:
        elemlist = nodes[i].where
        m.cons8.add(sum(m.x[elem,p] for p in m.PR-{0} for elem in elemlist) >= 2*m.y[i]) 
#----------------------------------------------------------        
    m.cons9 = ConstraintList()
    for i in m.LN:
        doflist = nodes[i].dof
        for dof in doflist:
            m.cons9.add(m.d[dof] <= m.y[i]*m.dmax)
            m.cons9.add(m.d[dof] >= -1*m.y[i]*m.dmax)
#----------------------------------------------------------
    m.cons10 = ConstraintList()        
    for i in m.LE:
        for j in range(i+1, len(m.LE)):
            seg1 = LineString([[elements[i].nodei.x, elements[i].nodei.y], [elements[i].nodej.x, elements[i].nodej.y]])
            seg2 = LineString([[elements[j].nodei.x, elements[j].nodei.y], [elements[j].nodej.x, elements[j].nodej.y]])
            int_pt = seg1.intersects(seg2)
            toucher = seg1.touches(seg2)
            if int_pt == True and toucher == False:
                for it in itertools.product(m.PR, m.PR):
                    m.cons10.add(m.x[i,it[0]] + m.x[j,it[1]] <= 1)
            else:
                Constraint.Skip
#    m.write('salamat','lp')
# %% Solving the MILP model
    if wheresol == 'PC':
        if sol == 'CPLEX':
            msolver = SolverFactory('CPLEX')
#            msolver.options['threads'] = 8
#            msolver.options['display'] = 2 #for each iteration
            msolver.options['timelimit'] = 36000
        elif sol == 'gurobi':
            msolver = SolverFactory('gurobi')
#            msolver.options['threads'] = 20
#            msolver.options['concurrentmip'] = 8
            msolver.options['timelim'] = 36000
            #msolver.options['iterlim'] = default no limit
        elif sol == 'XPRESS':
            msolver = SolverFactory('XPRESS')
#            msolver.options['miptol'] =1e-6  #integer feasibility tolerance; default = 5e-6]
            msolver.options['maxtime'] = 36000   # n < 0 ==> stop LP or MIP search after -n seconds,
                                             # n = 0 ==> no time limit (default),
                                             # n > 0 ==> for MIP problems, stop after n seconds if a feasible solution has been found;
                                                # otherwise continue until a feasible solution has been found.
#            msolver.options['threads'] = 8   # default number of threads to use: -1 = automatic choice (based on hardware),  n > 0 ==> use n threads
#            msolver.options['lpiterlimit'] = ... # simplex iteration limit; default = 2147483647 = 2^31 - 1
#            msolver.options['lplog'] = 100   # frequency of printing simplex iteration log; default = 100
        elif sol == 'octeract-engine':
            msolver = SolverFactory('octeract-engine')
            # msolver.options['LP_SOLVER'] = 'GUROBI'
            msolver.options['MILP_SOLVER'] = 'GUROBI'
        else:
            print('What else?')
                
        solution = msolver.solve(m, tee=True)
    else:
        solver_manager = SolverManagerFactory('neos')
        solution = solver_manager.solve(m, solver=sol)   
#    log_infeasible_constraints(m)
#    solution = solver_manager.solve(m, solver='cplex', options={'integrality':1e-09})
#----------------------------------------------------------
    for i in m.LE:
        for p in m.PR:
            if value(m.x[i,p]==1):
                dof1 = elements[i].dof
                dal = []
                for j in dof1:
                    dal.append(value(m.d[j]))
                KELT = np.dot(elements[i].transform.T,np.dot(elements[i].profile[p].lmat,elements[i].transform))
                m.Q[i] = np.dot(KELT,dal)
    weight = value(m.obj)
    print('\n ************** LP is done with status "{}"! Weight is "{}" and solver is "{}"**************.'.format(solution.solver.termination_condition, np.round(weight,4), sol))
#    print('-------------------------------------------( BEAMS )-------------------------------------------')
    X=[]
    for i in m.LE:
        for k in m.PR:
            if m.x[i,k]==1:
                X.append([i,k])
    toprint = pd.DataFrame(index = m.LE, columns = ['Profile','area','N1','V1','M1','N2','V2','M2'])
    for i in m.LE:
        for p in m.PR:
            if m.x[i,p]==1:
                force = value(m.Q[i])
                toprint.loc[i] = [p,elements[i].profile[p].area,force[0],force[1],force[2],force[3],force[4],force[5]]
    toprint.dropna(inplace=True)
    print(toprint)
#    print('-------------------------------------------( NODES )-------------------------------------------')
    toprint2 = pd.DataFrame(index = m.LN, columns = ['y','U1','U2','UR3'])
    for i in m.LN:
        dof = nodes[i].dof
        dis = np.round([m.d[dof[0]].value, m.d[dof[1]].value, m.d[dof[2]].value],5)
        toprint2.loc[i] = [m.y[i].value,dis[0],dis[1],dis[2]]
#    print(toprint2)
#    print('-------------------------------------------( REACTIONS )---------------------------------------')
    for nn in boundary:
        doffs = nodes[nn].dof
        print('RF in node {} is: in X direction: "{}", in Y direction: "{}" and rotation: "{}"'.format(nn,np.round(m.RF[doffs[0]].value,4),np.round(m.RF[doffs[1]].value,4),np.round(m.RF[doffs[2]].value,4)))
    # duration = 300; fres = 1000; winsound.Beep(fres, duration)    
    return(X, weight)