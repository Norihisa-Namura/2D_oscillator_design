from turtle import position
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import sys

"""
Description of functions:

- vector_field: Function for drawing and saving

- plot_vector_field: Function only for plotting

"""

def vector_field():
    # settings
    rc('text', usetex=True)
    plt.rcParams["font.family"] = "Times new roman"

    # load data
    X1 = np.loadtxt("data/X1.dat",delimiter=',')
    X2 = np.loadtxt("data/X2.dat",delimiter=',')
    F1_rep = np.loadtxt("data/F1_rep.dat",delimiter=',')
    F2_rep = np.loadtxt("data/F2_rep.dat",delimiter=',')
    F1_true = np.loadtxt("data/F1_true.dat",delimiter=',')
    F2_true = np.loadtxt("data/F2_true.dat",delimiter=',')
    p = np.loadtxt("data/p_rep.dat",delimiter=',')

    x_tick = np.loadtxt("data/x_tick.dat",delimiter=',')
    y_tick = np.loadtxt("data/y_tick.dat",delimiter=',')

    with open("data/class_name.txt",'r') as f:
        class_name = f.readline()
        class_name = class_name.replace('\n','')
            
    # draw vector field
    if class_name == "sync_cluster" or class_name == "star":
        fig = plt.figure(figsize=(6,5))
        plt.rcParams["font.size"] = 18
        plot_vector_field(fig,1,X1,X2,p,F1_rep,F2_rep,x_tick,y_tick)
    else:
        text = ["(a)","(b)"]
        fig = plt.figure(figsize=(6,10))
        plt.rcParams["font.size"] = 18
        clim,cticks = plot_vector_field(fig,2,X1,X2,p,F1_true,F2_true,x_tick,y_tick,text=text)
        plot_vector_field(fig,1,X1,X2,p,F1_rep,F2_rep,x_tick,y_tick,clim=clim,cticks=cticks,text=text)

    # save figure
    save_name = class_name + "_field"# + save_name
    plt.savefig("figs/" + save_name + ".pdf",bbox_inches='tight')

def plot_vector_field(fig,num,X1,X2,p,F1,F2,x_tick,y_tick,clim=None,cticks=None,text=None):
    C = np.log1p(np.sqrt(F1*F1 + F2*F2)) # maps to colormap 
    if clim == None:
        max_C = max(np.ravel(C))
        min_C = min(np.ravel(C))
        clim = [min_C,max_C]
        cticks = [min_C,(min_C+max_C)/2,max_C]

    scale = np.sqrt(F1*F1 + F2*F2) + 1e-10
    nU,nV = F1/scale,F2/scale

    if text == None:
        ax = fig.add_subplot(1,1,1)
    else:
        ax = plt.subplot2grid((30, 10), ((num-1)*16, 1), rowspan=12,colspan=8)

    ax.plot(p[0,:],p[1,:],'k')
    mappable = ax.quiver(X1, X2, nU, nV, C, cmap='rainbow', scale=40,width=0.005,headwidth=5,headlength=3,headaxislength=2,alpha=0.8,clim=clim)
    cbar = fig.colorbar(mappable=mappable,ax=ax,ticks=cticks)
    cbar.ax.set_yticklabels(['{:.4g}'.format(np.round(np.exp(cticks[0])-1,4)),'{:.4g}'.format(np.exp(cticks[1])-1),'{:.4g}'.format(np.exp(cticks[2])-1)])

    ax.set_xticks(x_tick)
    ax.set_yticks(y_tick)
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    if text != None:
        ax.text(0.5, -0.27, text[num-1], ha='center', transform=ax.transAxes)

    return clim,cticks

if __name__ == '__main__':
    vector_field()