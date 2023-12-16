from jpype import JPackage, startJVM, getDefaultJVMPath
import numpy
import sys
import numpy as np

import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt


def print_persent_of_run (destination_num,number_of_node):
    print((destination_num/number_of_node)*100,"%")
    pass

#readFloatsFile
def readFloatsFile(filename):
    "Read a 2D array of floats from a given file"
    with open(filename) as f:
        # Space separate numbers, one time step per line, each column is a variable
        array = []
        for line in f:
            # read all lines
            if (line.startswith("%") or line.startswith("#")):
                # Assume this is a comment line
                continue
            if (len(line.split()) == 0):
                # Line is empty
                continue
            array.append([float(x) for x in line.split()])
    return array

def read_data(file_address):
    # Our python data file readers are a bit of a hack, python users will do better on this:
    sys.path.append("~/Desktop/MyCodes/InformationTheory/JavaPackage/jlizier-jidt-64a7a80/demos/python")
    # Add JIDT jar library to the path
    jarLocation = "./infodynamics.jar"
    # Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
    startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + jarLocation)
    # 0. Load/prepare the data:
    dataRaw = readFloatsFile(file_address)
    return numpy.array(dataRaw)

def TE_Kraskov_matrix (data,from_node,to_node,time_column):
    # 1. Construct the calculator:
    calcClass = JPackage("infodynamics.measures.continuous.kraskov").TransferEntropyCalculatorKraskov
    calc = calcClass()
    # 2. Set any properties to non-default values:
    #calc.setProperty("NOISE_LEVEL_TO_ADD", "0.0")
    calc.setProperty("NOISE_LEVEL_TO_ADD", "1.0E-3")
    # 3. Initialise the calculator for (re-)use:
    calc.initialise()
    # 4. Supply the sample data:
    calc.setObservations(data[:,from_node+time_column],data[:,to_node+time_column])
    # 5. Compute the estimate:
    result = calc.computeAverageLocalOfObservations()
    return result

def calculate_matrix_information (number_of_node,source_num,name_file,time_column):
    file_address_input="./input_data/"+name_file+".txt"
    file_address_output="./output_data/"+name_file+".txt"
    data = read_data(file_address_input)
    file = open (file_address_output, 'w')
    file.write('From_1_to_A')
    for i in range(number_of_node):
        file.write('\t')
        file.write(str(i))
    file.write('\n')
    for from_node in range(number_of_node):
        file.write(str(from_node))
        for to_node in range(number_of_node):
            file.write('\t')
            if source_num==-1:
                result = TE_Kraskov_matrix (data,from_node,to_node,time_column)
                file.write("%.4f" %
                    (result))
                #print(str(from_node)+'---->'+str(to_node)+'='+str("%.4f" %result))
            else:
                if source_num==from_node or source_num==to_node:
                    result = TE_Kraskov_matrix (data,from_node,to_node,time_column)
                    file.write("%.4f" %
                        (result))
                    #print(str(from_node)+'---->'+str(to_node)+'='+str("%.4f" %result))
                else:
                    file.write("nan") 
        print_persent_of_run(from_node,number_of_node)
        file.write('\n')
    file.close()
    pass

def calculate_total_information (name_file,number_of_node):
    file_address_output="./output_data/"+name_file+".txt"
    data = np.loadtxt(file_address_output, skiprows=1, usecols=range(1, number_of_node+1))
    file_address_output_total_info="./output_data/"+name_file+"_total_info.txt"
    file_new = open (file_address_output_total_info, 'w')
    print("total information for each node")
    total=0
    for target_node in range(number_of_node):
        total_each_node=0
        for loop_all_conter in range(number_of_node):
            total_each_node+=data[target_node,loop_all_conter]
            total_each_node+=data[loop_all_conter,target_node]
        total_each_node-=data[target_node,target_node]# Because it is calculated twice
        file_new.write(str(target_node)+'\t'+str("%.4f" % total_each_node)+'\n')
        total+=total_each_node
        print(str(target_node)+" --> "+str("%.4f" % total_each_node))
    print("total information")
    print("I_total --> "+str("%.4f" % total))
    file_new.write('I_total\t'+str("%.4f" % total)+'\n')
    file_new.close()
    pass

def heatmap_plot(name_file,number_of_node):
    def heatmap(data, row_labels, col_labels, ax=None,
                cbar_kw=None, cbarlabel="", **kwargs):
        """
        Create a heatmap from a numpy array and two lists of labels.
        Parameters
        ----------
        data
            A 2D numpy array of shape (M, N).
        row_labels
            A list or array of length M with the labels for the rows.
        col_labels
            A list or array of length N with the labels for the columns.
        ax
            A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
            not provided, use current axes or create a new one.  Optional.
        cbar_kw
            A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
        cbarlabel
            The label for the colorbar.  Optional.
        **kwargs
            All other arguments are forwarded to `imshow`.
        """
        if ax is None:
            ax = plt.gca()
        if cbar_kw is None:
            cbar_kw = {}
        # Plot the heatmap
        im = ax.imshow(data, **kwargs)
        # Create colorbar
        cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
        cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
        # Show all ticks and label them with the respective list entries.
        ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
        ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)
        # Let the horizontal axes labeling appear on top.
        ax.tick_params(top=False, bottom=True,
                    labeltop=False, labelbottom=True)
        # Rotate the tick labels and set their alignment.
        '''plt.setp(ax.get_xticklabels(), rotation=-10, ha="right",
                rotation_mode="anchor")'''
        # Turn spines off and create white grid.
        ax.spines[:].set_visible(False)
        ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
        ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
        ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
        ax.tick_params(which="minor", bottom=False, left=False)

        return im, cbar
    def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                        textcolors=("black", "white"),
                        threshold=None, **textkw):
        """
        A function to annotate a heatmap.
        Parameters
        ----------
        im
            The AxesImage to be labeled.
        data
            Data used to annotate.  If None, the image's data is used.  Optional.
        valfmt
            The format of the annotations inside the heatmap.  This should either
            use the string format method, e.g. "$ {x:.2f}", or be a
            `matplotlib.ticker.Formatter`.  Optional.
        textcolors
            A pair of colors.  The first is used for values below a threshold,
            the second for those above.  Optional.
        threshold
            Value in data units according to which the colors from textcolors are
            applied.  If None (the default) uses the middle of the colormap as
            separation.  Optional.
        **kwargs
            All other arguments are forwarded to each call to `text` used to create
            the text labels.
        """
        if not isinstance(data, (list, np.ndarray)):
            data = im.get_array()
        # Normalize the threshold to the images color range.
        if threshold is not None:
            threshold = im.norm(threshold)
        else:
            threshold = im.norm(data.max())/2.
        # Set default alignment to center, but allow it to be
        # overwritten by textkw.
        kw = dict(horizontalalignment="center",
                verticalalignment="center")
        kw.update(textkw)
        # Get the formatter in case a string is supplied
        if isinstance(valfmt, str):
            valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)
        # Loop over the data and create a `Text` for each "pixel".
        # Change the text's color depending on the data.
        texts = []
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
                text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
                texts.append(text)
        return texts

    file_address_output="./output_data/"+name_file+".txt"
    data = np.loadtxt(file_address_output, skiprows=1, usecols=range(1, number_of_node+1))
    yـaxis = np.arange(number_of_node)
    xـaxis = np.arange(number_of_node)
    fig, ax = plt.subplots()
    im, cbar = heatmap(data, yـaxis, xـaxis, ax=ax,
                    cmap="YlGn", cbarlabel="Information")
    if number_of_node<=100:
        texts = annotate_heatmap(im, valfmt="{x:.3f}")
    plt.xlabel("To node", fontsize=10)
    plt.ylabel("From node", fontsize=10)
    fig.tight_layout()
    plt.subplots_adjust(top = 0.92, bottom=0.08,left=0.08)
    plt.gcf().set_size_inches(8, 6.5)# don't change it
    plt.savefig("./output_data/"+name_file+'.png',dpi=300)
    pass

def scatter_plot_for_source_from(name_file,source_num,number_of_node):
    file_address_output="./output_data/"+name_file+".txt"
    data = np.loadtxt(file_address_output, skiprows=1, usecols=range(1, number_of_node+1))
    #print(data[source_num])
    x=[i for i in range (number_of_node)]
    fig, ax = plt.subplots()
    plt.scatter(x,data[source_num])
    plt.xlabel("index of node", fontsize=10)
    plt.ylabel("information of node"+str(source_num), fontsize=10)
    plt.xlim(0, number_of_node)
    plt.ylim(0, 0.3)
    fig.tight_layout()
    plt.subplots_adjust(top = 0.92, bottom=0.08,left=0.1)
    plt.gcf().set_size_inches(8, 6.5)# don't change it
    plt.savefig("./output_data/from-node"+name_file+'-scatter.png',dpi=300)
    pass

def scatter_plot_for_source_to (name_file,source_num,number_of_node):
    file_address_output="./output_data/"+name_file+".txt"
    data = np.loadtxt(file_address_output, skiprows=1, usecols=range(1, number_of_node+1))
    #print(data[:,100])
    x=[i for i in range (number_of_node)]
    fig, ax = plt.subplots()
    plt.scatter(x,data[:,100])
    plt.xlabel("index of node", fontsize=10)
    plt.ylabel("information of node"+str(source_num), fontsize=10)
    plt.xlim(0, number_of_node)
    plt.ylim(0, 0.3)
    fig.tight_layout()
    plt.subplots_adjust(top = 0.92, bottom=0.08,left=0.1)
    plt.gcf().set_size_inches(8, 6.5)# don't change it
    plt.savefig("./output_data/to-node"+name_file+'-scatter.png',dpi=300)
    pass

def main():
    name_file="k=2.700000"#"two_sine_wave_with_shifted"
    time_column=1       # Hint1: If you have the time column, put the number 1, otherwise, put the number 0
    number_of_node=1000    # Hint2: The node index starts from 1 if first column is time
    source_num=100       # Hint3: if you want to calculate all node source_num=-1
    calculate_matrix_information (number_of_node,source_num,name_file,time_column)
    calculate_total_information (name_file,number_of_node) # Hint4: if you want to calculate total information
    #heatmap_plot (name_file,number_of_node)
    scatter_plot_for_source_from (name_file,source_num,number_of_node)
    scatter_plot_for_source_to (name_file,source_num,number_of_node)
    pass

if __name__=="__main__":
    main()