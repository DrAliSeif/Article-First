from jpype import JPackage, startJVM, getDefaultJVMPath
import numpy
import sys


def print_persent_of_run (destination_num,number_of_node):
    per_persent=number_of_node/100
    ##print("TE_Kraskov (KSG)(col_"+str(source_num)+" -> col_"+str(destination_num)+") = %.4f nats" %
        #(result))
    if destination_num%per_persent==0:
        print(destination_num/per_persent,"%")
    pass


def TE_Kraskov (source_arr,destination_arr,from_node,to_node,file):
    # 1. Construct the calculator:
    calcClass = JPackage("infodynamics.measures.continuous.kraskov").TransferEntropyCalculatorKraskov
    calc = calcClass()
    # 2. Set any properties to non-default values:
    #calc.setProperty("NOISE_LEVEL_TO_ADD", "0.0")
    #calc.setProperty("NOISE_LEVEL_TO_ADD", "1.0E-3")
    # 3. Initialise the calculator for (re-)use:
    calc.initialise()
    # 4. Supply the sample data:
    calc.setObservations(source_arr, destination_arr)
    # 5. Compute the estimate:
    result = calc.computeAverageLocalOfObservations()
    file.write(str(from_node)+'\t'+str(to_node)+'\t'+"%.4f" %
        (result))
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


def calculate_information (number_of_node,source_num,file_address_input,file_address_output):
    data = read_data(file_address_input)
    source_arr = data[:,source_num]
    file = open (file_address_output, 'w')    
    for destination_num in range (1,number_of_node+1,1):#1,1001,1
        destination_arr = data[:,destination_num]
        TE_Kraskov (source_arr,destination_arr,source_num,destination_num,file)
        file.write('\t')
        TE_Kraskov (destination_arr,source_arr,destination_num,source_num,file)
        file.write('\n')
        print_persent_of_run(destination_num,number_of_node)
    file.close()
    pass



def main():
    # Hint: The node index starts from 1
    number_of_node=1000
    source_num=100
    file_address_input="./input_data/new1k01mean10.txt"
    file_address_output="./output_data/noinfinitnew1k01mean10.txt"

    calculate_information(number_of_node,source_num,file_address_input,file_address_output)

    pass

if __name__=="__main__":
    main()