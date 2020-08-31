import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
# %% Drawing
def Draw_MILP(nodes, elements, X, W, TLP, solver):
## Drawing the Nodes ----------------------------
    fig, ax = plt.subplots()
    G = nx.Graph()
    pos = {}
    node_names = {}
    for i in nodes.keys():
        G.add_node(i)
        pos.update({i:[nodes[i].x,nodes[i].y]})
        node_names.update({i:nodes[i].name})
    nx.draw_networkx_nodes(G, pos, node_color = 'k', alpha = 1, node_size = 200, node_shape = 'o', linewidths = 0)
    nx.draw_networkx_labels(G, pos, node_names, font_size = 10, font_color='#FFFF00')
## Drawing the edges ----------------------------
    edge_labels = {}
    edge_widths = {}
    for i in X:
        i_pos1 = elements[i[0]].nodei.name
        i_pos2 = elements[i[0]].nodej.name
        G.add_edge(i_pos1, i_pos2)
        edge_labels.update({(i_pos1, i_pos2): elements[i[0]].name})
#            edge_labels.update({(i_pos1, i_pos2): str([np.round(barsdf.iloc[i]['force'],2),np.round(barsdf.iloc[i]['stress'],2), np.round(barsdf.iloc[i]['NLParea'],2)])})
        edge_widths.update({(i_pos1, i_pos2):4*elements[i[0]].profile[i[1]].area})
    edge_widths1 = list(edge_widths.values())
    nx.draw_networkx_edges(G, pos, edge_color = 'k', width = edge_widths1, ax = ax)
    nx.draw_networkx_edge_labels(G, pos, edge_labels = edge_labels, font_size = 10, font_color = 'r', alpha = 1, rotate = True, label_pos = 0.4)
    plt.axis('off')
    plt.suptitle('Weight:' + str(np.round(W,4)) + ' Stime:' + str(TLP) + ' Solver:' + solver, fontsize = 10)
    plt.show()

def Draw_MINLP(nodes, celements, Y, W, TNLP, solver):
## Drawing the Nodes ----------------------------
    fig, ax = plt.subplots()
    G = nx.Graph()
    pos = {}
    node_names = {}
    for i in nodes.keys():
        G.add_node(i)
        pos.update({i:[nodes[i].x,nodes[i].y]})
        node_names.update({i:nodes[i].name})
    nx.draw_networkx_nodes(G, pos, node_color = 'k', alpha = 1, node_size = 200, node_shape = 'o', linewidths = 0)
    nx.draw_networkx_labels(G, pos, node_names, font_size = 10, font_color='#FFFF00')
## Drawing the edges ----------------------------
    edge_labels1 = {}
    edge_widths = {}
    for i in celements.keys():
        if Y[i] != 0:
            i_pos1 = celements[i].nodei.name
            i_pos2 = celements[i].nodej.name
            G.add_edge(i_pos1, i_pos2)
            edge_labels1.update({(i_pos1, i_pos2): celements[i].name})
            edge_widths.update({(i_pos1, i_pos2):4*Y[i]})
    edge_widths1 = list(edge_widths.values())
    nx.draw_networkx_edges(G, pos, edge_color = 'k', width = edge_widths1, ax = ax)
    nx.draw_networkx_edge_labels(G, pos, edge_labels = edge_labels1, font_size = 10, font_color = 'r', alpha = 1, rotate = True, label_pos = 0.4)
    plt.axis('off')
    plt.suptitle('Weight:' + str(np.round(W,4)) + ' Stime:' + str(TNLP) + ' Solver:' + solver,fontsize = 10)
    plt.show()
def Draw_GS(nodes, elements):
## Drawing the Nodes ----------------------------
    fig, ax = plt.subplots()
    G = nx.Graph()
    pos = {}
    node_names = {}
#    node_shapes = {}
    node_colors = []
    for i in nodes.keys():
        G.add_node(i)
        pos.update({i:[nodes[i].x,nodes[i].y]})
        node_names.update({i:nodes[i].name})
        if nodes[i].tip == 1: #boundary
#            node_shapes.update({i: '^'})
            node_colors.append('b')
        elif nodes[i].tip == 2: #load point
#            node_shapes.update({i: '$\\uparrow$'})
            if nodes[i].load > 0:
                node_colors.append('g')
            else:
                node_colors.append('r')
        else:
#            node_shapes.update({i: 'o'})
            node_colors.append('k')
#    print(node_shapes)   
    nx.draw_networkx_nodes(G, pos, node_color = node_colors , alpha = 1, node_size = 100, node_shape = 'o', linewidths = 0)
#    nx.draw_networkx_labels(G, pos, node_names, font_size = 10, font_color='#FFFF00')
## Drawing the edges ----------------------------
    edge_labels1 = {}
    edge_widths = {}
    for i in elements.keys():
        i_pos1 = elements[i].nodei.name
        i_pos2 = elements[i].nodej.name
        G.add_edge(i_pos1, i_pos2)
        edge_labels1.update({(i_pos1, i_pos2): elements[i].name})
        edge_widths.update({(i_pos1, i_pos2):1})
    edge_widths1 = list(edge_widths.values())
    nx.draw_networkx_edges(G, pos, edge_color = 'k', width = edge_widths1, ax = ax)
#    nx.draw_networkx_edge_labels(G, pos, edge_labels = edge_labels1, font_size = 10, font_color = 'r', alpha = 1, rotate = True, label_pos = 0.4)
    plt.axis('off')
    plt.show()