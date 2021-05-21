import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import numpy as np
from shapely.geometry import LineString
from pyomo.environ import Param, ConcreteModel, Var, Objective, ConstraintList, value, minimize, Binary, Constraint
from pyomo.opt import SolverFactory, SolverStatus, TerminationCondition
import itertools
import time

def warm(E, nodes, elements, r_set, dmax, smax):
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
            load_node_dof.append(nodes[i].dof[0])
            load_node_dof.append(nodes[i].dof[1])
            load_value.append(nodes[i].load)
    f = np.ravel([nodes[i].f for i in nodes.keys()])

    m = ConcreteModel()
    m.LN = set(nodes.keys())  # List of nodes
    m.LE = set(elements.keys())  # List of elements
    m.PR = set(range(len(r_set)))  # List of profiles
    m.dofs = set(np.ravel([nodes[i].dof for i in nodes.keys()]))
    m.nfree = nfree
    m.free = m.dofs - m.nfree
    m.rho = 1
    m.E = E
    m.dmax = dmax
    m.Smax = smax

    m.d = Var(m.dofs, initialize=0)  # displacement along each degrees of freedom
    for i in m.nfree:  # fixed nodes have all displacements equal to zeros
        m.d[i].fix(0)
    m.x = Var(m.LE, m.PR, domain=Binary, initialize=0)  # index i is used for the elements or beams
    m.y = Var(m.LN, domain=Binary, initialize=0)  # indicator of each nodes to be connected to the structure
    m.v = Var(m.LE, {0, 1, 2}, m.PR,
              initialize=0)  # Elongation or contraction of beam i with profile p in three dimensions
    m.s = Var(m.LE, {0, 1, 2}, initialize=0)  # forces (axial, shear, rotation) on each element
    m.Q = Param(m.LE, initialize=0, mutable=True)  # forces (axial or rotational) on each degree of freedom
    m.RF = Param(m.nfree, initialize=0,
                 mutable=True)  # Reaction forces on the boundary condition fixed degrees of freedom
    m.M = 100000  # Big M for the purpose of linearization of logic constraint
    # --------------------------------------------------------------OBJECTIVE---------------------------------------------------------------------------
    m.obj = Objective(expr=sum([elements[i].profile[p].vol * m.x[i, p] for i in m.LE for p in m.PR]), sense=minimize)
    # -------------------------------------------------------------CONSTRAINTS---------------------------------------------------------------------------
    m.cons1 = ConstraintList()  # Equilibrium equations
    for satr in m.free:
        temp1 = 0
        for i in m.LE:
            for p in m.PR:
                temp1 += (elements[i].profile[p].KE[0]) * (m.v[i, 0, p]) * (elements[i].B[0][satr]) + \
                         (elements[i].profile[p].KE[1]) * (m.v[i, 1, p]) * (elements[i].B[1][satr]) + \
                         (elements[i].profile[p].KE[2]) * (m.v[i, 2, p]) * (elements[i].B[2][satr])
        m.cons1.add(temp1 == f[satr])
    # ----------------------------------------------------------
    m.cons2 = ConstraintList()  # When there is a beam it's elongation should be equal to the difference between displacement in its two nodes.
    for i in m.LE:
        for j in {0, 1, 2}:
            m.cons2.add(
                sum(m.v[i, j, p] for p in m.PR - {0}) - sum(elements[i].B[j][d] * m.d[d] for d in m.dofs) <= m.M * (
                        1 - sum(m.x[i, p] for p in m.PR - {0})))
            m.cons2.add(sum(m.v[i, j, p] for p in m.PR - {0}) - sum(
                elements[i].B[j][d] * m.d[d] for d in m.dofs) >= -1 * m.M * (1 - sum(m.x[i, p] for p in m.PR - {0})))
    # ----------------------------------------------------------
    m.cons3 = ConstraintList()  # Force on each beam should be equal to elongation * KE (scalar-3 values for each profile)
    for i in m.LE:
        for j in {0, 1, 2}:
            m.cons3.add(elements[i].profile[p].KE[j] * (sum(m.v[i, j, p] for p in m.PR - {0})) - m.s[i, j] <= m.M * (
                    1 - sum(m.x[i, p] for p in m.PR - {0})))
            m.cons3.add(
                elements[i].profile[p].KE[j] * (sum(m.v[i, j, p] for p in m.PR - {0})) - m.s[i, j] >= -1 * m.M * (
                        1 - sum(m.x[i, p] for p in m.PR - {0})))
    # ----------------------------------------------------------
    m.cons4 = ConstraintList()  # Constraint to correlate cross-sectional area and bearable force
    for i in m.LE:
        for j in {0, 1, 2}:
            m.cons4.add(m.s[i, j] <= m.M * sum(m.x[i, p] for p in m.PR - {0}))
            m.cons4.add(m.s[i, j] >= -m.M * sum(m.x[i, p] for p in m.PR - {0}))
    # ----------------------------------------------------------
    m.cons5 = ConstraintList()  # Correlate existence of a beam with the force on that
    for i in m.LE:
        for j in {0, 1, 2}:
            m.cons5.add(m.v[i, j, p] <= m.M * sum(m.x[i, p] for p in m.PR - {0}))
            m.cons5.add(m.v[i, j, p] >= -m.M * sum(m.x[i, p] for p in m.PR - {0}))
    # ----------------------------------------------------------
    m.cons6 = ConstraintList()  # Each beam can only have one cross sectional area from the candidate list
    for i in m.LE:
        m.cons6.add(sum(m.x[i, p] for p in m.PR - {0}) <= 1)
    # ----------------------------------------------------------
    m.cons7 = ConstraintList()  # If x binary for a beam equal to 1, y binary for both end nodes should be equal to 1.
    for i in m.LE:
        m.cons7.add(
            sum(m.x[i, p] for p in m.PR - {0}) <= 0.5 * (m.y[elements[i].nodei.name] + m.y[elements[i].nodej.name]))
    # ----------------------------------------------------------
    m.cons8 = ConstraintList()  # Degree of each node (except load and boundary nodes) should be greater than 1 to remove the hanging beams.
    for i in m.LN - boundary - load_node:
        elemlist = nodes[i].where
        m.cons8.add(sum(m.x[elem, p] for p in m.PR - {0} for elem in elemlist) >= 2 * m.y[i])
    # ----------------------------------------------------------
    m.cons9 = ConstraintList()  # Absoulute value of displacement in each DoF should be less than max possible displacement
    for i in m.LN - boundary:
        doflist = nodes[i].dof
        for dof in doflist:
            m.cons9.add(m.d[dof] <= m.y[i] * m.dmax)
            m.cons9.add(m.d[dof] >= -m.y[i] * m.dmax)
    # ----------------------------------------------------------
    m.cons10 = ConstraintList()  # Removing over-crossing in the beams.
    for i in m.LE:
        for j in range(i + 1, len(m.LE)):
            seg1 = LineString([[elements[i].nodei.x, elements[i].nodei.y], [elements[i].nodej.x, elements[i].nodej.y]])
            seg2 = LineString([[elements[j].nodei.x, elements[j].nodei.y], [elements[j].nodej.x, elements[j].nodej.y]])
            int_pt = seg1.intersects(seg2)
            toucher = seg1.touches(seg2)
            if int_pt == True and toucher == False:
                for it in itertools.product(m.PR, m.PR):
                    m.cons10.add(m.x[i, it[0]] + m.x[j, it[1]] <= 1)
            else:
                Constraint.Skip
    msolver = SolverFactory('gurobi_persistent')
    msolver.options['TimeLimit'] = 300  # Time limit increase from 10800 to 18000 on 11.23.2020
    msolver.options['LogToConsole'] = 1
    msolver.options['DisplayInterval'] = 10
    msolver.options['Heuristics'] = 1
    msolver.options['MIPFocus'] = 1
    msolver.options['FeasibilityTol'] = 1e-3
    msolver.options['MIPGap'] = 0.01
    msolver.options['Threads'] = 16
    msolver.options['IntegralityFocus'] = 1
    msolver.options['BranchDir'] = -1
    msolver.options['Cuts'] = 3
    msolver.set_instance(m)
    Wstart = time.time()
    solution = msolver.solve(m, tee=True, keepfiles=True)#
    TWS = np.round(time.time() - Wstart, 3)
    m.solutions.store_to(solution)
    if (solution.solver.status == SolverStatus.ok) and (
            solution.solver.termination_condition == TerminationCondition.optimal):
        print('WS is optimal and solution pickled')
    elif solution.solver.termination_condition == TerminationCondition.infeasible:
        print('OOPS the WS is infeasible')
    else:
        print('***Another try with more time on WS, time limit increased to 1200 sec.***')
        msolver.options['TimeLimit'] = 900
        msolver.options['FeasibilityTol'] = 1e-2
        solution = msolver.solve(m, tee=True, keepfiles=True)
        print('Solver Status: {}'.format(solution.solver.status))
    weight = value(m.obj)
    X = {}
    Y = {}
    D = {}
    V = {}
    S = {}
    for i in m.LE:
        for k in m.PR:
            X.update({(i, k): m.x[i, k].value})
    for i in m.LE:
        for j in {0, 1, 2}:
            S.update({(i, j): m.s[i, j].value})
    for i in m.LE:
        for j in {0, 1, 2}:
            for p in m.PR:
                V.update({(i, j, p): m.v[i, j, p].value})
    for i in m.LN:
        Y.update({i: m.y[i].value})
    for i in m.dofs:
        D.update({i: m.d[i].value})
    return X, Y, D, V, S, weight, TWS
