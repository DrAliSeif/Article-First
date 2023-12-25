import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn

def read_data(filename_address):
    file_address_output="./data_file_input/"+filename_address+".txt"
    data = np.loadtxt(file_address_output)
    #print(pd.DataFrame(data))
    return data,len(data[:,0]),len(data[0,:])

def plot_se_heatmap(name_file,matrix):
    if len(matrix)<=10:
        ax = sn.heatmap(matrix, vmin=-1, vmax=1,annot=True, fmt=".1f", linewidth=.5, cmap="coolwarm")#,annot=True, fmt=".1f" ---> for show number amount
    else:
        ax = sn.heatmap(matrix, vmin=-1, vmax=1 ,cmap="seismic")#,annot=True, fmt=".1f" ---> for show number amount
    #ax.set(xlabel="", ylabel="")
    ax.xaxis.tick_top()
    plt.title("node", fontsize=10)
    plt.ylabel("node", fontsize=10)
    #plt.clabel("Correlation", fontsize=10)
    plt.xlim([0, len(matrix)])
    plt.ylim([len(matrix),0])
    plt.setp(ax.get_xticklabels(), rotation=90, ha="left",rotation_mode="anchor")
    plt.subplots_adjust(top = 0.90, bottom=0.02,left=0.1)
    plt.gcf().set_size_inches(8, 6.5)# don't change it
    plt.savefig("./output_data/"+name_file+".png",dpi=300)
    pass

def correlation(one_snapshot_of_data,number_of_node):
    matrix=np.zeros((number_of_node, number_of_node))
    for i in range (number_of_node):
        for j in range (number_of_node):
            matrix[i,j]=np.cos(one_snapshot_of_data[j]-one_snapshot_of_data[i])
    return matrix

def main():
    name_file="k=0.000000"
    data,rows,column=read_data(name_file)
    matrix_total=np.zeros((column, column))
    #rows=100
    for target in range (rows):
        print((target+1)/rows)
        matrix_total+=correlation(data[target,:],column)
    matrix_total=matrix_total/rows
    plot_se_heatmap(name_file,matrix_total)
    pass

if __name__=="__main__":
    main()