""" Created on Tue May 12 19:02:12 2020 
    Ground structure for the discrete model
    OOP is used in this version of the code
    @author: ozt0008 
    Modified on Jun 18 to cover both discrete and continuous models
    Modified on Nov 25"""

import numpy as np

dim = 2  # dimension of the problem
ndof = 3  # Number of degrees of freedom for each node
# distance = 10 # Square grid length(distance of each of the nodes from its neighbour nodes)

__all__ = ['Generate', 'Node', 'Element', 'Profile', 'CElement']


# --------------------------------------------------------------------------------------
class Generate:  # Ground structure generator
    #    r_set = [] # set of available discrete circular cross-sectional radius
    def __init__(self, M, N, fix_nodes, load_nodes, load_values, E, r_set):
        global distance
        distance = 50 / (M - 1)
        self.n_nodes = M * N
        self.nodes = {}  # Dictionary of Node instances to be filled
        self.elements = {}  # Dictionary of discrete Element instances to be filled
        self.celements = {}  # Dictionary of continuous Element instances to be filled
        for i in range(self.n_nodes):
            if i in load_nodes:
                l_n = 2
                l_v = load_values[load_nodes.index(i)]
            elif i in fix_nodes:
                l_n = 1
                l_v = 0
            else:
                l_n = 0
                l_v = 0
            self.nodes[i] = Node(i, i % M, i // M, l_n, l_v)
        elem = 0
        for i in range(self.n_nodes):
            for j in range(i + 1, self.n_nodes):
                temp_element = Element(self.nodes[i], self.nodes[j], M, N, elem, E, r_set)
                temp_celement = CElement(self.nodes[i], self.nodes[j], M, N, elem, E)
                if temp_element.okay == 1:
                    self.nodes[i].where.append(elem)
                    self.nodes[j].where.append(elem)
                    self.elements[elem] = temp_element
                    #                    print(temp_element.profile[0].KE)
                    self.celements[elem] = temp_celement
                    #                    print(temp_celement.KE)
                    elem += 1


# --------------------------------------------------------------------------------------
class Node:
    def __init__(self, node_name, x, y, tip, load):
        self.name = node_name
        self.x = x * distance
        self.y = y * distance
        self.tip = tip
        self.load = load
        self.dof = self._dof
        self.where = []
        self.f = self._f

    @property
    def _dof(self):
        s = self.name * ndof
        return [s, s + 1, s + 2]

    @property
    def _f(self):
        if self.tip != 2:
            return [0, 0, 0]
        else:
            return [0, self.load, 0]


# --------------------------------------------------------------------------------------
class Element:
    def __init__(self, nodei, nodej, M, N, elem, E, r_set):
        self.nodei = nodei
        self.nodej = nodej
        self.orient = [self.nodei.name, self.nodej.name]
        self.dof = np.concatenate([nodei.dof, nodej.dof], axis=None)
        self.length = self._length
        self.theta = self._theta
        self.cosan = self._cosan
        self.sinan = self._sinan
        self.transform = np.array([[self.cosan, self.sinan, 0, 0, 0, 0],
                                   [-self.sinan, self.cosan, 0, 0, 0, 0],
                                   [0, 0, 1, 0, 0, 0],
                                   [0, 0, 0, self.cosan, self.sinan, 0],
                                   [0, 0, 0, -self.sinan, self.cosan, 0],
                                   [0, 0, 0, 0, 0, 1]])
        self.Tmat = np.zeros([M * N * ndof, 6])
        for i in range(6):
            self.Tmat[self.dof[i], i] = 1

        if ((nodei.x == nodej.x and np.abs(nodei.y - nodej.y) > distance) or (np.abs(nodei.x - nodej.x) > distance)):
            self.okay = -1
        else:
            self.okay = 1
            self.name = elem
            self.profile = {}
            for i in range(len(r_set)):
                self.profile[i] = Profile(self, E, r_set[i])
            self.B = self.profile[1].B

    @property
    def _length(self):
        return np.sqrt((self.nodej.y - self.nodei.y) ** 2 + (self.nodej.x - self.nodei.x) ** 2)

    @property
    def _theta(self):
        return np.arctan2((self.nodej.y - self.nodei.y), (self.nodej.x - self.nodei.x))

    @property
    def _cosan(self):
        return (self.nodej.x - self.nodei.x) / self.length

    @property
    def _sinan(self):
        return (self.nodej.y - self.nodei.y) / self.length


# --------------------------------------------------------------------------------------
class Profile:
    def __init__(self, el, E, r):
        self.el = el
        if r != 0:
            self.r = r
            self.area = np.pi * (r ** 2)
            self.I = (self.area ** 2) / (4 * np.pi)
            self.b1g = np.dot(el.Tmat, np.dot(el.transform.T, np.array([[-1, 0, 0, 1, 0, 0]]).T))
            self.b2g = np.dot(el.Tmat,
                              np.dot(el.transform.T, np.array([[0, 2 / el.length, 1, 0, -2 / el.length, 1]]).T))
            self.b3g = np.dot(el.Tmat, np.dot(el.transform.T, np.array([[0, 0, -1, 0, 0, 1]]).T))
            self.ke1 = (E * self.area) / el.length
            self.ke2 = (3 * E * self.I) / el.length
            self.ke3 = (E * self.I) / el.length
            self.KE = [self.ke1, self.ke2, self.ke3]
            self.vol = self.area * el.length
            self.smat = self._smatrix
            self.B = self._B
            self.lmat = self._lmatrix
        else:
            self.r = r
            self.area = 0
            self.I = 0
            self.b1g = np.dot(el.Tmat, np.dot(el.transform.T, np.array([[0, 0, 0, 0, 0, 0]]).T))
            self.b2g = np.dot(el.Tmat, np.dot(el.transform.T, np.array([[0, 0, 0, 0, 0, 0]]).T))
            self.b3g = np.dot(el.Tmat, np.dot(el.transform.T, np.array([[0, 0, 0, 0, 0, 0]]).T))
            self.ke1 = 0
            self.ke2 = 0
            self.ke3 = 0
            self.KE = [self.ke1, self.ke2, self.ke3]
            self.vol = 0
            self.smat = self._smatrix
            self.B = self._B
            self.lmat = self._lmatrix

    @property
    def _smatrix(self):
        temp = np.zeros([self.b1g.size, self.b1g.size])
        dim1 = self.b1g.size
        b1 = [i * self.ke1 for i in np.dot(self.b1g, self.b1g.T)]
        b2 = [i * self.ke2 for i in np.dot(self.b2g, self.b2g.T)]
        b3 = [i * self.ke3 for i in np.dot(self.b3g, self.b3g.T)]
        for i in range(dim1):
            for j in range(dim1):
                temp[i][j] += b1[i][j] + b2[i][j] + b3[i][j]
        return temp

    @property
    def _lmatrix(self):
        temp1 = np.zeros([6, 6])
        if self.r != 0:
            b1loc = np.dot(np.array([[-1, 0, 0, 1, 0, 0]]).T, np.array([[-1, 0, 0, 1, 0, 0]]))
            b2loc = np.dot(np.array([[0, 2 / self.el.length, 1, 0, -2 / self.el.length, 1]]).T,
                           np.array([[0, 2 / self.el.length, 1, 0, -2 / self.el.length, 1]]))
            b3loc = np.dot(np.array([[0, 0, -1, 0, 0, 1]]).T, np.array([[0, 0, -1, 0, 0, 1]]))
            for i in range(6):
                for j in range(6):
                    temp1[i][j] += self.ke1 * b1loc[i][j] + self.ke2 * b2loc[i][j] + self.ke3 * b3loc[i][j]
        return temp1

    @property
    def _B(self):
        B = [self.b1g.T.flatten(), self.b2g.T.flatten(), self.b3g.T.flatten()]
        return B
    # --------------------------------------------------------------------------------------


class CElement:
    def __init__(self, nodei, nodej, M, N, elem, E):
        self.nodei = nodei
        self.nodej = nodej
        self.orient = [self.nodei.name, self.nodej.name]
        self.dof = np.concatenate([nodei.dof, nodej.dof], axis=None)
        self.length = self._length
        self.theta = self._theta
        self.cosan = self._cosan
        self.sinan = self._sinan
        self.transform = np.array([[self.cosan, self.sinan, 0, 0, 0, 0],
                                   [-self.sinan, self.cosan, 0, 0, 0, 0],
                                   [0, 0, 1, 0, 0, 0],
                                   [0, 0, 0, self.cosan, self.sinan, 0],
                                   [0, 0, 0, -self.sinan, self.cosan, 0],
                                   [0, 0, 0, 0, 0, 1]])
        self.Tmat = np.zeros([M * N * ndof, 6])
        for i in range(6):
            self.Tmat[self.dof[i], i] = 1
        self.name = elem
        self.b1g = np.dot(self.Tmat, np.dot(self.transform.T, np.array([[-1, 0, 0, 1, 0, 0]]).T))
        self.b2g = np.dot(self.Tmat,
                          np.dot(self.transform.T, np.array([[0, 2 / self.length, 1, 0, -2 / self.length, 1]]).T))
        self.b3g = np.dot(self.Tmat, np.dot(self.transform.T, np.array([[0, 0, -1, 0, 0, 1]]).T))
        self.ke1 = E / self.length
        self.ke2 = (3 * E) / (4 * np.pi * self.length)
        self.ke3 = E / (4 * np.pi * self.length)
        self.KE = [self.ke1, self.ke2, self.ke3]
        #        self.smat = self._smatrix
        self.B = self._B
        self.bloc = self._Bloc

    #    @property
    #    def _smatrix(self):
    #        ctemp = np.zeros([self.b1g.size,self.b1g.size])
    #        dim1 = self.b1g.size
    #        b1 = [i*self.ke1 for i in np.dot(self.b1g, self.b1g.T)]
    #        b2 = [i*self.ke2 for i in np.dot(self.b2g, self.b2g.T)]
    #        b3 = [i*self.ke3 for i in np.dot(self.b3g, self.b3g.T)]
    #        for i in range(dim1):
    #            for j in range(dim1):
    #                ctemp[i][j] += b1[i][j] + b2[i][j] + b3[i][j]
    #        return ctemp
    @property
    def _B(self):
        B = [self.b1g.T.flatten(), self.b2g.T.flatten(), self.b3g.T.flatten()]
        return B

    @property
    def _Bloc(self):
        b1loc = np.dot(np.array([[-1, 0, 0, 1, 0, 0]]).T, np.array([[-1, 0, 0, 1, 0, 0]]))
        b2loc = np.dot(np.array([[0, 2 / self.length, 1, 0, -2 / self.length, 1]]).T,
                       np.array([[0, 2 / self.length, 1, 0, -2 / self.length, 1]]))
        b3loc = np.dot(np.array([[0, 0, -1, 0, 0, 1]]).T, np.array([[0, 0, -1, 0, 0, 1]]))
        return [b1loc, b2loc, b3loc]

    @property
    def _length(self):
        return np.sqrt((self.nodej.y - self.nodei.y) ** 2 + (self.nodej.x - self.nodei.x) ** 2)

    @property
    def _theta(self):
        return np.arctan2((self.nodej.y - self.nodei.y), (self.nodej.x - self.nodei.x))

    @property
    def _cosan(self):
        return (self.nodej.x - self.nodei.x) / self.length

    @property
    def _sinan(self):
        return (self.nodej.y - self.nodei.y) / self.length
# print(*dir(model...), sep = '\n')
