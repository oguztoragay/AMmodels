###### Updated on 09/12/2020 for DoE
###### Updated on 01/27/2021 for Paper
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import os

def Draw_Warm(nodes, elements, X, TWS, W, load, fname, foldername):
    Xkey = []
    for i in elements.keys():
        if X[i, 1] == 1:
            Xkey.append(i)
    Xval = [0.5 for i in Xkey]
    Xdic = dict(zip(Xkey, Xval))
    nodeset1 = [elements[i].nodei.name for i in Xkey]
    nodeset2 = [elements[i].nodej.name for i in Xkey]
    node_list = list(set(nodeset1) | set(nodeset2))
    ## Drawing the Nodes ----------------------------
    fig, ax = plt.subplots(figsize=(5, 5))
    G = nx.Graph()
    pos = {}
    node_names = {}
    node_colors = []
    for i in nodes.keys():
        G.add_node(i)
        pos.update({i: [nodes[i].x, nodes[i].y]})
        node_names.update({i: nodes[i].name})
        if nodes[i].tip == 1:  # boundary
            node_colors.append('b')
        elif nodes[i].tip == 2:  # load point
            if nodes[i].load > 0:
                node_colors.append('g')
            else:
                node_colors.append('r')
        else:
            if i in node_list:
                node_colors.append('k')
            else:
                node_colors.append('lightgrey')
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, alpha=1, node_size=200, node_shape='o', linewidths=0)
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
            edge_widths.update({(i_pos1, i_pos2): 10})
            edge_colors.update({(i_pos1, i_pos2): 'r'})
            edge_styles.update({(i_pos1, i_pos2): 'solid'})
        else:
            edge_widths.update({(i_pos1, i_pos2): 1})
            edge_colors.update({(i_pos1, i_pos2): 'lightgrey'})
            edge_styles.update({(i_pos1, i_pos2): 'dashed'})
    edge_colors1 = list(edge_colors.values())
    edge_styles1 = list(edge_styles.values())
    edge_widths1 = list(edge_widths.values())
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors1, width=edge_widths1, ax=ax, style=edge_styles1)
    # nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=10, font_color='r', alpha=1,
    #                              rotate=True,
    #                              label_pos=0.4)
    plt.axis('off')
    # plt.suptitle('|Stime:' + str(TWS) + '|Weight:' + str(np.round(W, 4)) + '|Load:' + str(load), fontsize=11)
    plt.show()
    try:
        os.mkdir(foldername)
    except OSError as error:
        print(error)
    fig.savefig(str(foldername+'/'+fname+'.pdf'), bbox_inches='tight', pad_inches=0.01)

def Draw_MILP(nodes, elements, X, S, W, TLP, ii, jj, solver, dmax, fname, foldername):
    Xkey = [i[0] for i in X]
    Xval = [i[1] for i in X]
    Sval = [np.round(i, 1) for i in S]
    Xdic = dict(zip(Xkey, Xval))
    Sdic = dict(zip(Xkey, Sval))
    nodeset1 = [elements[i].nodei.name for i in Xkey]
    nodeset2 = [elements[i].nodej.name for i in Xkey]
    node_list = list(set(nodeset1) | set(nodeset2))
    ## Drawing the Nodes ----------------------------
    fig, ax = plt.subplots(figsize=(5, 5))
    G = nx.Graph()
    pos = {}
    node_names = {}
    node_colors = []
    for i in nodes.keys():
        G.add_node(i)
        pos.update({i: [nodes[i].x, nodes[i].y]})
        node_names.update({i: nodes[i].name})
        if nodes[i].tip == 1:  # boundary
            node_colors.append('b')
        elif nodes[i].tip == 2:  # load point
            if nodes[i].load > 0:
                node_colors.append('g')
            else:
                node_colors.append('r')
        else:
            if i in node_list:
                node_colors.append('k')
            else:
                node_colors.append('lightgrey')
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, alpha=1, node_size=200, node_shape='o', linewidths=0)
    # nx.draw_networkx_labels(G, pos, node_names, font_size=10, font_color='#FFFF00')
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
            edge_labels.update({(i_pos1, i_pos2): str(elements[i].profile[Xdic[i]].r)})
            edge_widths.update({(i_pos1, i_pos2): 40 * elements[i].profile[Xdic[i]].r})
            edge_colors.update({(i_pos1, i_pos2): 'k'})
            edge_styles.update({(i_pos1, i_pos2): 'solid'})
        else:
            edge_widths.update({(i_pos1, i_pos2): 1})
            edge_colors.update({(i_pos1, i_pos2): 'lightgrey'})
            edge_styles.update({(i_pos1, i_pos2): 'dashed'})
    edge_colors1 = list(edge_colors.values())
    edge_styles1 = list(edge_styles.values())
    edge_widths1 = list(edge_widths.values())
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors1, width=edge_widths1, ax=ax, style=edge_styles1)
    # nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=10, font_color='r', alpha=1, rotate=True,
    #                              label_pos=0.4)
    plt.axis('off')
    # plt.suptitle('|Dmax:' + str(dmax) + '|Solver:' + solver + '|Stime:' + str(TLP) + '|Weight:' + str(
    #     np.round(W, 4)) + '\nCS:' + str(jj) + '|Load:' + str(ii), fontsize=11)
#    plt.suptitle('Model: ' + str('CS'+str(a+1)) +
#                 '\nTime: ' + str(TLP) +
#                 '\nWeight: ' + str(np.round(W, 2)), fontsize=12, ha='left', x=0.1)
    plt.show()
#    parent_dir='C:/Users/ozt0008/Documents/OneDrive - Auburn University/1 AM/Models/model PYTHON/7 January 2021/01.27.2021/'
#    path = os.path.join(parent_dir, foldername)
#    try:
#        os.mkdir(foldername)
#    except OSError as error:
#        print(error)
    fig.savefig(str(foldername+'/'+fname+'.pdf'), bbox_inches='tight', pad_inches=0)

def Draw_MINLP(nodes, celements, Y, W, TNLP, ii, solver, dmax, fname, foldername):
    Ydic = {i: Y[i] for i in range(0, len(Y))}
    nodeset1 = []
    nodeset2 = []
    for i in Ydic:
        if Ydic[i] >= 0.1256:
            nodeset1.append(celements[i].nodei.name)
            nodeset2.append(celements[i].nodej.name)
    node_list = list(set(nodeset1) | set(nodeset2))
    ## Drawing the Nodes ----------------------------
    fig, ax = plt.subplots(figsize=(5, 5))
    G = nx.Graph()
    pos = {}
    node_names = {}
    node_colors = []
    for i in nodes.keys():
        G.add_node(i)
        pos.update({i: [nodes[i].x, nodes[i].y]})
        node_names.update({i: nodes[i].name})
        if nodes[i].tip == 1:  # boundary
            node_colors.append('b')
        elif nodes[i].tip == 2:  # load point
            if nodes[i].load > 0:
                node_colors.append('g')
            else:
                node_colors.append('r')
        else:
            if i in node_list:
                node_colors.append('k')
            else:
                node_colors.append('lightgrey')
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, alpha=1, node_size=200, node_shape='o', linewidths=0)
    # nx.draw_networkx_labels(G, pos, node_names, font_size=10, font_color='#FFFF00')
    ## Drawing the edges ----------------------------
    edge_labels = {}
    edge_widths = {}
    edge_colors = {}
    edge_styles = {}
    for i in celements.keys():
        i_pos1 = celements[i].nodei.name
        i_pos2 = celements[i].nodej.name
        G.add_edge(i_pos1, i_pos2)
        if Ydic[i] >= 0.1256:
            edge_labels.update({(i_pos1, i_pos2): str(np.round(np.sqrt(Y[i] / np.pi), 3))})#[celements[i].name,
            edge_widths.update({(i_pos1, i_pos2): 40 * (np.sqrt(Y[i] / np.pi))})
            edge_colors.update({(i_pos1, i_pos2): 'k'})
            edge_styles.update({(i_pos1, i_pos2): 'solid'})
        else:
            edge_widths.update({(i_pos1, i_pos2): 1})
            edge_colors.update({(i_pos1, i_pos2): 'lightgrey'})
            edge_styles.update({(i_pos1, i_pos2): 'dashed'})
    edge_colors1 = list(edge_colors.values())
    edge_styles1 = list(edge_styles.values())
    edge_widths1 = list(edge_widths.values())
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors1, width=edge_widths1, ax=ax, style=edge_styles1)
    # nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=10, font_color='r', alpha=1, rotate=True, label_pos=0.4)
    plt.axis('off')
    # plt.suptitle('|Dmax:' + str(dmax) + '|Solver:' + solver + '|Stime:' + str(TNLP) + '|Weight:' + str(
    #     np.round(W, 4)) + '\nLoad:' + str(ii), fontsize=11)
#    plt.suptitle('Model: ' + solver +
#                 '\nTime: ' + str(TNLP) +
#                 '\nWeight: ' + str(np.round(W, 2)), fontsize=12, ha='left', x=0.1)
    plt.show()
    fig.savefig(str(foldername+'/'+fname+'.pdf'), bbox_inches='tight', pad_inches=0)

def Draw_GROUND_solid(nodes, elements, fname, foldername):
    fig, ax = plt.subplots(figsize=(5, 5))
    G = nx.Graph()
    pos = {}
    node_names = {}
    for i in nodes.keys():
        G.add_node(i)
        pos.update({i: [nodes[i].x, nodes[i].y]})
        node_names.update({i: nodes[i].name})
    nx.draw_networkx_nodes(G, pos, node_color='k', alpha=1, node_size=200, node_shape='o', linewidths=0)
    nx.draw_networkx_labels(G, pos, node_names, font_size=10, font_color='#FFFF00')
    for i in elements.keys():
        i_pos1 = elements[i].nodei.name
        i_pos2 = elements[i].nodej.name
        G.add_edge(i_pos1, i_pos2)
    nx.draw_networkx_edges(G, pos, edge_color='k', width=2, ax=ax, style='solid')
    plt.axis('off')
    plt.show()
    fig.savefig(str(foldername + '/' + fname + '.pdf'), bbox_inches='tight', pad_inches=0.01)

def Draw_GROUND_dashed(nodes, elements, fname, foldername):
    fig, ax = plt.subplots(figsize=(5, 5))
    G = nx.Graph()
    colors = []
    pos = {}
    for i in nodes.keys():
        pos.update({i: [nodes[i].x, nodes[i].y]})
        if nodes[i].tip == 1:  # boundary
            G.add_node(i, s='o')
            colors.append('b')
        elif nodes[i].tip == 2:  # load point
            G.add_node(i, s='o')
            if nodes[i].load > 0:
                x_load = nodes[i].x
                y_load = nodes[i].y
                ax.arrow(x_load, y_load, 0, 6, width=0.7, facecolor='k', head_width=2, length_includes_head=False)
                colors.append('g')
            else:
                x_load = nodes[i].x
                y_load = nodes[i].y
                ax.arrow(x_load, y_load, 0, -6, width=0.7, facecolor='k', head_width=2, length_includes_head=False)
                colors.append('r')
        else:
            G.add_node(i, s='o')
            colors.append('k')
    nodeShapes = set((aShape[1]["s"] for aShape in G.nodes(data=True)))
    for aShape in nodeShapes:
        nx.draw_networkx_nodes(G, pos, node_shape=aShape, node_color=colors,
                               nodelist=[sNode[0] for sNode in filter(lambda x: x[1]["s"] == aShape, G.nodes(data=True))],
                               alpha=1, node_size=150,  linewidths=0, )#node_color='k'
    for i in elements.keys():
        i_pos1 = elements[i].nodei.name
        i_pos2 = elements[i].nodej.name
        G.add_edge(i_pos1, i_pos2)
    nx.draw_networkx_edges(G, pos, edge_color='lightgrey', width=1, ax=ax, style='dashed')
    plt.axis('off')
    plt.show()
    fig.savefig(str(foldername + '/' + fname + '.pdf'), bbox_inches='tight', pad_inches=0)
