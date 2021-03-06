import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import numpy as np
from shapely.geometry import LineString
from pyomo.environ import Param, ConcreteModel, Var, Objective, ConstraintList, value, minimize, Binary, Constraint
from pyomo.opt import SolverFactory
from pyomo.opt.parallel import SolverManagerFactory
from pyomo.util.infeasible import log_infeasible_constraints
import pandas as pd
import time

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# %% Sets, Variables and Parameters
def NLPpyo(E, nodes, celements, r2_set, dmax, smax, sol, wheresol):
    boundary = set()
    load_node = set()
    nfree = set()
    load_node_dof = []
    load_value = []
    for i in nodes.keys():
        if nodes[i].tip == 1:
            boundary.add(i)
            nfree.update(nodes[i].dof)
        elif nodes[i].tip == 2:
            load_node.add(i)
            load_node_dof.append(nodes[i].dof[1])
            load_value.append(nodes[i].load)
    f = np.ravel([nodes[i].f for i in nodes.keys()])
    ## --------------------------------------------------------------------------------------------------------------------------------------------------
    m = ConcreteModel()
    m.LN = set(nodes.keys())
    m.LE = set(celements.keys())
    m.dofs = set(np.ravel([nodes[i].dof for i in nodes.keys()]))
    m.nfree = nfree
    m.free = m.dofs - m.nfree
    m.rho = 1
    m.E = E
    m.dmax = dmax
    m.amin = np.pi * (r2_set[0] ** 2)
    m.amax = np.pi * (r2_set[2] ** 2)
    m.Smax = smax
    ## --------------------------------------------------------------------------------------------------------------------------------------------------
    def lerule(m, i):
        return celements[i].length

    m.d = Var(m.dofs, initialize=0)
    for i in m.nfree:
        m.d[i].fix(0)
    m.a = Var(m.LE, initialize=r2_set[1])
    m.x = Var(m.LE, domain=Binary, initialize=1)
    m.y = Var(m.LN, domain=Binary, initialize=0)
    m.v = Var(m.LE, {0, 1, 2})  # Elongation or contraction of beam i , initialize=0
    m.RF = Var(m.nfree, initialize=0)  # Reactions on the boundary nodes
    m.le = Param(m.LE, rule=lerule)
    m.M = Param(initialize=10000)

    # %% Generating Objective function and Constraints
    def obj_rule(m):
        return sum(m.a[i] * m.le[i] for i in m.LE)

    m.z = Objective(rule=obj_rule, sense=minimize)
    # -------------------------------------------------------------CONSTRAINTS---------------------------------------------------------------------------
    m.cons1 = ConstraintList()
    for satr in m.free:
        m.ps = set(nodes[np.floor(satr / 3)].where)
        temp1 = 0
        for i in m.LE:
            temp1 += (celements[i].KE[0] * m.a[i]) * (
                    celements[i].B[0][satr] * sum(celements[i].B[0][di0] * m.d[di0] for di0 in m.dofs)) + \
                     (celements[i].KE[1] * m.a[i] ** 2) * (
                             celements[i].B[1][satr] * sum(celements[i].B[1][di0] * m.d[di0] for di0 in m.dofs)) + \
                     (celements[i].KE[2] * m.a[i] ** 2) * (
                             celements[i].B[2][satr] * sum(celements[i].B[2][di0] * m.d[di0] for di0 in m.dofs))
        # if satr in m.nfree:
        #     m.cons1.add(temp1 == m.RF[satr])
        # else:
        m.cons1.add(temp1 == f[satr])
    # -----------------------------------------------------------------------------------------
    m.cons2 = ConstraintList()
    for i in m.LE:
        m.cons2.add(m.a[i] <= m.amax * m.x[i])
        m.cons2.add(m.a[i] >= m.amin * m.x[i])
    # -----------------------------------------------------------------------------------------
    m.cons3 = ConstraintList()
    for i in m.LN - boundary - load_node:
        elemlist = nodes[i].where
        m.cons3.add(2 * m.y[i] <= (sum(m.x[elem] for elem in elemlist)))
    # -----------------------------------------------------------------------------------------
    m.cons4 = ConstraintList()
    for i in m.free:
        m.cons4.add(m.d[i] <= m.dmax * m.y[i // 3])
        m.cons4.add(m.d[i] >= -m.dmax * m.y[i // 3])
    # -----------------------------------------------------------------------------------------
    m.cons5 = ConstraintList()
    for i in m.LE:
        m.cons5.add(2 * m.x[i] <= m.y[celements[i].nodei.name] + m.y[celements[i].nodej.name])
    # -----------------------------------------------------------------------------------------
    m.cons6 = ConstraintList()
    for i in m.LE:
        for j in range(i + 1, len(m.LE)):
            seg1 = LineString(
                [[celements[i].nodei.x, celements[i].nodei.y], [celements[i].nodej.x, celements[i].nodej.y]])
            seg2 = LineString(
                [[celements[j].nodei.x, celements[j].nodei.y], [celements[j].nodej.x, celements[j].nodej.y]])
            int_pt = seg1.intersects(seg2)
            toucher = seg1.touches(seg2)
            if int_pt == True and toucher == False:
                m.cons6.add(m.x[i] + m.x[j] <= 1)
            else:
                Constraint.Skip
    # -----------------------------------------------------------------------------------------
    # %% Solving MINLP model
	solver = SolverFactory('BARON')
    solver.options['MaxIter'] = -1
    solver.options['threads'] = 16
    solver.options['MaxTime'] = 7200
    solver.options['numsol'] = 1
    solver.options['outlev'] = 1
    solver.options['prtime'] = 100
    solver.options['lpsolver'] = 'cplex'
    solver.options['nlpsol'] = 4
    NLPstart = time.time()
    solution = solver.solve(m, tee=True)
    TNLP = np.round(time.time()-NLPstart, 3)
    data1_minlp = solution.Problem._list
    data2_minlp = solution.solver._list
    weight1 = value(m.z)
    print('\n ************** NLP is done with status "{}"! Weight is "{}" and solver is "{}"**************.'.format(
            solution.solver.termination_condition, np.round(weight1, 4), sol))
    # -----------------------------------------------------------------------------------------
    force_loc = {}
    for i in m.LE:
        dofz = celements[i].dof
        d1 = m.d[dofz[0]].value
        d2 = m.d[dofz[1]].value
        d3 = m.d[dofz[2]].value
        d4 = m.d[dofz[3]].value
        d5 = m.d[dofz[4]].value
        d6 = m.d[dofz[5]].value
        dz = [d1, d2, d3, d4, d5, d6]
        KEL = np.multiply(celements[i].ke1 * (value(m.a[i])), celements[i].bloc[0]) + \
              np.multiply(celements[i].ke2 * (value(m.a[i])) ** 2, celements[i].bloc[1]) + \
              np.multiply(celements[i].ke3 * (value(m.a[i])) ** 2, celements[i].bloc[2])
        KELT = np.dot(celements[i].transform.T, np.dot(KEL, celements[i].transform))
        F = np.dot(KELT, dz)
        force_loc[i] = F.T
    print('-------------------------------------------( BEAMS )-------------------------------------------')
    Y = []
    for i in m.LE:
        Y.append(m.a[i].value)
    toprint = pd.DataFrame(index=m.LE, columns=['No', 'Area', 'Radius', 'X', 'N1', 'V1', 'M1', 'N2', 'V2', 'M2'])
    for i in m.LE:
        toprint.loc[i] = [i, value(m.a[i]), np.sqrt(value(m.a[i]) / np.pi), value(m.x[i]), force_loc[i][0],
                          force_loc[i][1], force_loc[i][2], force_loc[i][3], force_loc[i][4], force_loc[i][5]]
    indexNames = toprint[(toprint['Area'] == 0)].index
    toprint.drop(indexNames, inplace=True)
    print(toprint)
    print('-------------------------------------------( NODES )-------------------------------------------')
    toprint2 = pd.DataFrame(index=nodes.keys(), columns=['y', 'U1', 'U2', 'UR3'])
    for i in nodes.keys():
        dofz = nodes[i].dof
        dis = [m.d[dofz[0]].value, m.d[dofz[1]].value, m.d[dofz[2]].value]
        toprint2.loc[i] = [m.y[i].value, dis[0], dis[1], dis[2]]
    print(toprint2)
    print('-------------------------------------------( REACTIONS )---------------------------------------')
    for nn in boundary:
        doffs = nodes[nn].dof
        print('RF in node {} is: in X direction: "{}", in Y direction: "{}" and rotation: "{}"'.format(nn, np.round(
            m.RF[doffs[0]].value, 4), np.round(m.RF[doffs[1]].value, 4), np.round(m.RF[doffs[2]].value, 4)))
    return Y, weight1, TNLP, data1_minlp, data2_minlp
