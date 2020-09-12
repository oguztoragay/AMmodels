###### Updated on 09/12/2020 for DoE
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
# %% Drawing
def Draw_MILP(nodes, elements, X, W, TLP, ii, jj, solver,dmax):
    Xkey = [i[0] for i in X]
    Xval = [i[1] for i in X]
    Xdic = dict(zip(Xkey,Xval))
    nodeset1 = [elements[i].nodei.name for i in Xkey]
    nodeset2 = [elements[i].nodej.name for i in Xkey]
    node_list = list(set(nodeset1) | set(nodeset2))
## Drawing the Nodes ----------------------------
    fig, ax = plt.subplots()
    G = nx.Graph()
    pos = {}
    node_names = {}
    node_colors = []
    for i in nodes.keys():
        G.add_node(i)
        pos.update({i:[nodes[i].x,nodes[i].y]})
        node_names.update({i:nodes[i].name})
        if nodes[i].tip == 1: #boundary
            node_colors.append('b')
        elif nodes[i].tip == 2: #load point
            if nodes[i].load > 0:
                node_colors.append('g')
            else:
                node_colors.append('r')
        else:
            if i in node_list:
                node_colors.append('k')
            else:
                node_colors.append('gray')
    nx.draw_networkx_nodes(G, pos, node_color = node_colors, alpha = 1, node_size = 200, node_shape = 'o', linewidths = 0)
    nx.draw_networkx_labels(G, pos, node_names, font_size = 10, font_color='#FFFF00')
## Drawing the edges ----------------------------
    edge_labels = {}
    edge_widths = {}
    edge_colors = {}
    edge_styles = {}
    for i in elements.keys():
        i_pos1 = elements[i].nodei.name
        i_pos2 = elements[i].nodej.name
        G.add_edge(i_pos1, i_pos2)
        if i in Xdic.keys():
            edge_labels.update({(i_pos1, i_pos2): str([i, Xdic[i]])})
            edge_widths.update({(i_pos1, i_pos2):10*elements[i].profile[Xdic[i]].area})
            edge_colors.update({(i_pos1, i_pos2):'k'})
            edge_styles.update({(i_pos1, i_pos2):'solid'})
        else:
            edge_widths.update({(i_pos1, i_pos2):1})
            edge_colors.update({(i_pos1, i_pos2):'gray'})
            edge_styles.update({(i_pos1, i_pos2):'dashed'})
    edge_colors1 = list(edge_colors.values())
    edge_styles1 = list(edge_styles.values())
    edge_widths1 = list(edge_widths.values())
    nx.draw_networkx_edges(G, pos, edge_color = edge_colors1, width = edge_widths1, ax = ax, style = edge_styles1)
    nx.draw_networkx_edge_labels(G, pos, edge_labels = edge_labels, font_size = 10, font_color = 'r', alpha = 1, rotate = True, label_pos = 0.4)
    plt.axis('off')
    plt.suptitle('|Dmax:' + str(dmax) + '|Solver:' + solver + '|Stime:' + str(TLP) + '|Weight:' + str(np.round(W,4)) + '\nCS:' + str(jj) + '|Load:' + str(ii), fontsize = 11)
    plt.show()

#def Draw_GS(nodes, elements):
### Drawing the Nodes ----------------------------
#    fig, ax = plt.subplots()
#    G = nx.Graph()
#    pos = {}
#    node_names = {}
#    node_colors = []
#    for i in nodes.keys():
#        G.add_node(i)
#        pos.update({i:[nodes[i].x,nodes[i].y]})
#        node_names.update({i:nodes[i].name})
#        if nodes[i].tip == 1: #boundary
#            node_colors.append('b')
#        elif nodes[i].tip == 2: #load point
#            if nodes[i].load > 0:
#                node_colors.append('g')
#            else:
#                node_colors.append('r')
#        else:
#            node_colors.append('k')  
#    nx.draw_networkx_nodes(G, pos, node_color = node_colors , alpha = 1, node_size = 100, node_shape = 'o', linewidths = 0)
### Drawing the edges ----------------------------
#    edge_labels1 = {}
#    edge_widths = {}
#    for i in elements.keys():
#        i_pos1 = elements[i].nodei.name
#        i_pos2 = elements[i].nodej.name
#        G.add_edge(i_pos1, i_pos2)
#        edge_labels1.update({(i_pos1, i_pos2): elements[i].name})
#        edge_widths.update({(i_pos1, i_pos2):1})
#    edge_widths1 = list(edge_widths.values())
#    nx.draw_networkx_edges(G, pos, edge_color = 'k', width = edge_widths1, ax = ax)
#    plt.axis('off')
#    plt.show()
    
def Draw_MINLP(nodes, celements, Y, W, TNLP, solver):
    Ydic = {i: Y[i] for i in range(0, len(Y))}
    nodeset1 = []
    nodeset2 = []
    for i in Ydic:
        if Ydic[i] !=0:
            nodeset1.append(celements[i].nodei.name)
            nodeset2.append(celements[i].nodej.name)
    node_list = list(set(nodeset1) | set(nodeset2))
## Drawing the Nodes ----------------------------
    fig, ax = plt.subplots()
    G = nx.Graph()
    pos = {}
    node_names = {}
    node_colors = []
    for i in nodes.keys():
        G.add_node(i)
        pos.update({i:[nodes[i].x,nodes[i].y]})
        node_names.update({i:nodes[i].name})
        if nodes[i].tip == 1: #boundary
            node_colors.append('b')
        elif nodes[i].tip == 2: #load point
            if nodes[i].load > 0:
                node_colors.append('g')
            else:
                node_colors.append('r')
        else:
            if i in node_list:
                node_colors.append('k')
            else:
                node_colors.append('gray')
    nx.draw_networkx_nodes(G, pos, node_color = node_colors, alpha = 1, node_size = 200, node_shape = 'o', linewidths = 0)
    nx.draw_networkx_labels(G, pos, node_names, font_size = 10, font_color='#FFFF00')
## Drawing the edges ----------------------------
    edge_labels = {}
    edge_widths = {}
    edge_colors = {}
    edge_styles = {}
    for i in celements.keys():
        i_pos1 = celements[i].nodei.name
        i_pos2 = celements[i].nodej.name
        G.add_edge(i_pos1, i_pos2)
        if Ydic[i]!=0:
            edge_labels.update({(i_pos1, i_pos2): celements[i].name})
            edge_widths.update({(i_pos1, i_pos2):10*Y[i]})
            edge_colors.update({(i_pos1, i_pos2):'k'})
            edge_styles.update({(i_pos1, i_pos2):'solid'})
        else:
            edge_widths.update({(i_pos1, i_pos2):1})
            edge_colors.update({(i_pos1, i_pos2):'gray'})
            edge_styles.update({(i_pos1, i_pos2):'dashed'})
    edge_colors1 = list(edge_colors.values())
    edge_styles1 = list(edge_styles.values())
    edge_widths1 = list(edge_widths.values())
    nx.draw_networkx_edges(G, pos, edge_color = edge_colors1, width = edge_widths1, ax = ax, style = edge_styles1)
    nx.draw_networkx_edge_labels(G, pos, edge_labels = edge_labels, font_size = 10, font_color = 'r', alpha = 1, rotate = True, label_pos = 0.4)
    plt.axis('off')
    plt.suptitle('Weight:' + str(np.round(W,4)) + ' Stime:' + str(TNLP) + ' Solver:' + solver,fontsize = 10)
    plt.show()