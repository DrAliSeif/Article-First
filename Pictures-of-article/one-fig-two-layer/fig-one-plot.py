import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
import matplotlib.patheffects as path_effects
from mpl_toolkits.mplot3d.art3d import Line3DCollection


def main():
    cols = ['#3333ff', '#e60029']
    number_of_node=10
    # Imagine you have three node-aligned snapshots of a network
    G1 = nx.random_geometric_graph(number_of_node, 1, seed=896803)
    G2 = nx.random_geometric_graph(number_of_node, 1, seed=896803)
    #G3 = nx.karate_club_graph()

    pos = nx.circular_layout(nx.random_geometric_graph(number_of_node, 1, seed=896803)) # assuming common node location
    graphs = [G1,G2]#, G3]
    w = 18
    h = 18

    fig, ax = plt.subplots(1, 1, figsize=(w,h), dpi=100, subplot_kw={'projection':'3d'})

    for gi, G in enumerate(graphs):
        # node positions
        xs = list(list(zip(*list(pos.values())))[0])
        ys = list(list(zip(*list(pos.values())))[1])
        zs = [gi]*len(xs) # set a common z-position of the nodes 

        # node colors
        cs = [cols[gi]]*len(xs)
        
        # if you want to have between-layer connections
        if gi > 0:
            thru_nodes = np.arange(number_of_node)
            lines3d_between = [(list(pos[i])+[gi-1],list(pos[i])+[gi]) for i in thru_nodes]
            between_lines = Line3DCollection(lines3d_between, zorder=gi, color='.1',
                                            alpha=0.4, linestyle='--', linewidth=3)
            ax.add_collection3d(between_lines)

        # add within-layer edges 
        lines3d = [(list(pos[i])+[gi],list(pos[j])+[gi]) for i,j in G.edges()]
        line_collection = Line3DCollection(lines3d, zorder=gi-1, color=cols[gi], alpha=1, linewidth=3)
        ax.add_collection3d(line_collection)
        
        # now add nodes
        ax.scatter(xs, ys, zs,color=cs, s=990,  marker='o', alpha=1, zorder=gi+1, linewidth=3, edgecolor='k')
        
        # if you want labels...
        '''for li, lab in enumerate(list(G.nodes())):
            ax.text(xs[li], ys[li], zs[li]+0.007, lab, color='#262626', zorder=gi+200, fontsize=37,font='Calibri',
                    ha='center', va='center')'''
        
        # add a plane to designate the layer
        xdiff = max(xs)-min(xs)
        ydiff = max(ys)-min(ys)
        ymin = min(ys)-ydiff*0.1
        ymax = max(ys)+ydiff*0.1
        xmin = min(xs)-xdiff*0.1 * (w/h)
        xmax = max(xs)+xdiff*0.1 * (w/h)
        xx, yy = np.meshgrid([xmin, xmax],[ymin, ymax])
        zz = np.zeros(xx.shape)+gi
        ax.plot_surface(xx, yy, zz, color=cols[gi], alpha=0.1, zorder=gi)

        # add label layer
        if gi==0:
            stingit="1"
        else:
            stingit="2"
        layertext = ax.text(2.2, 2.8, gi*0.95+0.5, " layer"+stingit,
                            color=cols[gi], fontsize=45, zorder=1e5, ha='left', va='center')
        

    # set them all at the same x,y,zlims
    ax.set_ylim(min(ys)-ydiff*0.1,max(ys)+ydiff*0.1)
    ax.set_xlim(min(xs)-xdiff*0.1,max(xs)+xdiff*0.1)
    ax.set_zlim(-0.1, len(graphs) - 1 + 0.1)

    # select viewing angle
    angle = 30
    height_angle = 200
    ax.view_init(height_angle, angle)

    # how much do you want to zoom into the fig
    ax.dist = 9.0

    ax.set_axis_off()
    plt.subplots_adjust(top = 0.97, bottom=0.08, hspace=0.24, wspace=0.34)
    plt.gcf().set_size_inches(20, 17)
    #plt.savefig("fig-v1.png", dpi=100)
    plt.savefig("fig-v1.eps")

    pass

if __name__=="__main__":
    main()