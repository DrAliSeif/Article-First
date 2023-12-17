import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

def Read_Data(address):return np.loadtxt(address)

def plot_framework(number_of_pic):
    ax = plt.subplot(2,2,number_of_pic)
    plt.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=True, labeltop=False, pad=15)
    plt.tick_params(axis='y', which='both' , left=True, right=True, labelbottom=True, labeltop=False, pad=15)
    plt.xticks(fontsize=31,color= '#262626')
    plt.yticks(fontsize=31,color= '#262626')
    return ax


def plot_ran (ax,data):
    ax.plot(data[:,0], data[:,1] ,'o-', color="#3333ff",linewidth=2, markersize=14, markerfacecolor='none')
    ax.plot(data[:,0], data[:,2] , linestyle='dashed', marker='o', color="#3333ff",linewidth=2, markersize=8)
    ax.set_ylim(0,1)
    ax.set_xlim(0,3)
    ax.set_xlabel('Intralayer coupling strength ($\sigma$)',  fontsize=31)
    ax.set_ylabel('Synchronization ($\mathrm{r}^\mathrm{a}$)',  fontsize=31)
    ax.tick_params(axis='both', which='both', bottom=True, top=True, left=True, right=True)
    ax.tick_params(axis='both', direction='in', which='major', length=22)
    ax.tick_params(axis='both', direction='in', which='minor', length=6)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    pass


def plot_insert(ax,data,final,start,final_y):
    final100=int(final*100)
    start100=int(start*100)
    X = data[start100:final100, 0]
    Y = data[start100:final100, 1]
    Y2 = data[start100:final100, 2]
    # insert
    ax.plot(X, Y,'o-', color="#3333ff",linewidth=1, markersize=8, markerfacecolor='none')
    ax.plot(X, Y2, linestyle='dashed', marker='o', color="#3333ff",linewidth=1, markersize=4)
    ax.set_ylim(0,final_y)
    ax.set_xlim(start,final)
    ax.tick_params(axis='both', which='both', bottom=True, top=True, left=True, right=True)
    ax.tick_params(axis='both', direction='in', which='major', length=12)
    ax.tick_params(axis='both', direction='in', which='minor', length=4)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    pass

def plot_framework2():
    plt.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=True, labeltop=False, pad=15)
    plt.tick_params(axis='y', which='both' , left=True, right=True, labelbottom=True, labeltop=False, pad=15)
    plt.xticks(fontsize=20,color= '#262626')
    plt.yticks(fontsize=20,color= '#262626')
    pass

def text_plot():
    plt.text(0.97,1.58, '(a)',  fontsize=30,
         bbox=dict(boxstyle="round",
                   ec="#262626",
                   fc=(1, 1, 1, 0.80),
                   ))
    plt.text(3.51,1.58, '(b)',  fontsize=30,
         bbox=dict(boxstyle="round",
                   ec="#262626",
                   fc=(1, 1, 1, 0.80),
                   ))
    plt.text(0.97,-0.84, '(c)',  fontsize=30,
         bbox=dict(boxstyle="round",
                   ec="#262626",
                   fc=(1, 1, 1, 0.80),
                   ))
    plt.text(3.51,-0.84, '(d)',  fontsize=30,
         bbox=dict(boxstyle="round",
                   ec="#262626",
                   fc=(1, 1, 1, 0.80),
                   ))
    pass

def main():
    L_U=Read_Data('./a/r0_0to3.txt')
    R_U=Read_Data('./b/r70_0to3.txt')
    L_D=Read_Data('./c/s0_0to3.txt')
    R_D=Read_Data('./d/s70fb.txt')
    fig = plt.figure()

    axes=[plot_framework(i) for i in range(1, 5)]
    plot_ran(axes[0],L_U)
    plot_ran(axes[1],R_U)
    plot_ran(axes[2],L_D)
    plot_ran(axes[3],R_D)


    axes2 = fig.add_axes([0.23, 0.13, 0.18, 0.2]) # inset axes
    plot_framework2()
    plot_insert(axes2,L_D,0.5,0,1)

    axes3 = fig.add_axes([0.23, 0.614, 0.18, 0.2]) # inset axes
    plot_framework2()
    plot_insert(axes3,L_U,1,0,1)


    '''axes4 = fig.add_axes([0.73, 0.13, 0.18, 0.2]) # inset axes
    plot_framework2()
    plot_insert(axes4,R_D,0.8,0.3,0.5)
    '''
    text_plot()

    plt.subplots_adjust(top = 0.97, bottom=0.08, hspace=0.2, wspace=0.44)
    plt.gcf().set_size_inches(22, 17)
    plt.savefig("fig3_100_new.png", dpi=100)
    pass

if __name__=="__main__":
    main()