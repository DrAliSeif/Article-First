from jpype import JPackage, startJVM, getDefaultJVMPath
import numpy
import sys
import numpy as np

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
    #calc.setProperty("NOISE_LEVEL_TO_ADD", "1.0E-3")
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

def calculate_total_information (name_file,number_of_node,source_num):
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

def main():
    name_file="two_sine_wave_with_shifted_noise"#"two_sine_wave_with_shifted"
    time_column=0       # Hint1: If you have the time column, put the number 1, otherwise, put the number 0
    number_of_node=2    # Hint2: The node index starts from 1 if first column is time
    source_num=-1       # Hint3: if you want to calculate all node source_num=-1
    calculate_matrix_information (number_of_node,source_num,name_file,time_column)
    calculate_total_information (name_file,number_of_node,source_num) # Hint4: if you want to calculate total information
    pass

if __name__=="__main__":
    main()