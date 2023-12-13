#include"Kuramoto.Version3.h"//import library Kuramoto
#include <time.h>

int main() {//Beginning main
//data.txt and Example file
//Read and definition data,Number_of_node,Phases_initial,frequency,adj,copling,delay,time
    const double* data=read_data(new double[count_rows_file("data.txt")]);//read data from data.txt and write them to pointer 1d
    const int Number_of_node = int(double(data[1]));//Number_of_node=N=1000
    const double* frequency_layer1 = read_initial_1D("W=Natural frequency/Layer1", Number_of_node);//frequency  W
    const double* frequency_layer2 = read_initial_1D("W=Natural frequency/Saika0.8a", Number_of_node);
    const int* const* adj_layer1 = read_initial_2D("A=Intralayer adjacency matrix/Layer1", Number_of_node);//adjacency matrix  A
    const int* const* adj_layer2 = read_initial_2D("A=Intralayer adjacency matrix/Layer2", Number_of_node);
    const double* Phases_initial_layer1 = read_initial_1D("P=Initial Phases/Phases_initial_layer2_origin", Number_of_node);//Initial Phases  P
    const double* Phases_initial_layer2 = read_initial_1D("P=Initial Phases/Phases_initial_layer2_origin", Number_of_node);
    cout << "|----------------------------------------------------------------------------------------------------|\n" << endl;
    
    //definitions
    //const double Degree_Radian=1.567;
    const double delay = double(data[8]);
    const int history = int(delay / data[3]) + 1;//Create history for delay arrey
    const double landa=10;
    const int time_stationary = int(data[4] * 0.2);//example T=20 time_stationary= 10
    const int Number_Steps_time_stationary = int(time_stationary / data[3]);//for example T=20 dt=0.01 >> = 1000

    //pragma omp parallel for
    // initial
    //omp_set_num_threads(20);
    //#pragma omp parallel for
    for (int Degree = 1570; Degree <= 1570; Degree+=3) {
        //double frequency=double(double(frequency_10)/10.0);
        double Degree_Radian=Degree/1000.;
        ostringstream ostrDegree_Radian;// declaring output string stream
        ostrDegree_Radian << Degree_Radian;// Sending a number as a stream output
        string strDegree_Radian = ostrDegree_Radian.str();// the str() converts number into string
        ofstream Average_Syncrony("Save/Average_Syncrony(couplig_SyncL1_SyncL2)/Degree_Radian="+strDegree_Radian+".txt");//Create Syncrony file
        double copling = double(data[5]);
        double* Phases_next_layer1 = new double[Number_of_node];// Definition Phases next
        double* Phases_next_layer2 = new double[Number_of_node];
        double** Phases_history_delay_layer1 = Hystorydelay_phases(Number_of_node, history, Phases_initial_layer1);// history
        double** Phases_history_delay_layer2 = Hystorydelay_phases(Number_of_node, history, Phases_initial_layer2);// delay
        double* Phases_layer1_previous = previous_phases(Number_of_node, history, Phases_history_delay_layer1);//Phases changer
        double* Phases_layer2_previous = previous_phases(Number_of_node, history, Phases_history_delay_layer2);//[node][delay]
        // copling loop
        while (copling < (data[7])) {
            time_t start = time(NULL);
            ostringstream ostrcopling;// declaring output string stream
            ostrcopling << copling;// Sending a number as a stream output
            string strcopling = ostrcopling.str();// the str() converts number into string
            ofstream Phases_layer2("Save/Phases/Degree_Radian="+strDegree_Radian+"_copling="+strcopling+"layer2(time)VS(Node).txt");//     ---
            ofstream Phases_layer1("Save/Phases/Degree_Radian="+strDegree_Radian+"_copling="+strcopling+"layer1(time)VS(Node).txt");
            double Total_syncrony_layer1 = 0;
            double Total_syncrony_layer2 = 0;
            // time loop
            double time_loop = double(data[2]);// reset time for new time
            while (time_loop < (data[4] + data[3])) {
                // Runge Kutta and scaling phase
                Connected_Constant_Runge_Kutta_4(Degree_Radian, Number_of_node, data[3], delay, copling, frequency_layer1, adj_layer1, Phases_layer1_previous, Phases_history_delay_layer1, Phases_layer2_previous,Phases_next_layer1,landa);
                Connected_Constant_Runge_Kutta_4(Degree_Radian, Number_of_node, data[3], delay, copling, frequency_layer2, adj_layer2, Phases_layer2_previous, Phases_history_delay_layer2, Phases_layer1_previous,Phases_next_layer2,landa);
                scale_pi(Number_of_node, Phases_layer1_previous);// change values to be in range 0 to 2*Pi
                scale_pi(Number_of_node, Phases_layer2_previous);
                double syncrony_layer1 = order_parameter(Number_of_node, Phases_layer1_previous);// order parameters
                double syncrony_layer2 = order_parameter(Number_of_node, Phases_layer2_previous);// order parameters
                Phases_layer2 << time_loop << '\t';//Phases_layer1 << time << '\t';
                Phases_layer1 << time_loop << '\t';
                for (int i = 0; i < Number_of_node; i++) {
                    Phases_layer2 << Phases_layer2_previous[i] << '\t';//Phases_layer1 << Phases_layer1_previous[i] << '\t';
                    Phases_layer1 << Phases_layer1_previous[i] << '\t';
                }
                Phases_layer2 << endl;// Phases_layer1 << endl;
                Phases_layer1 << endl;
                if (time_loop > time_stationary) {// add sync to total sync
                    Total_syncrony_layer1 += syncrony_layer1;
                    Total_syncrony_layer2 += syncrony_layer2;
                }
                time_loop += data[3];
            }
            Total_syncrony_layer1=Total_syncrony_layer1/(Number_Steps_time_stationary*4);// calculate total sync and pint it
            Total_syncrony_layer2=Total_syncrony_layer2/(Number_Steps_time_stationary*4);
            time_t end = time(NULL);
            cout<< copling << '\t' << Total_syncrony_layer2 <<'\t' <<"Execution Time: "<< (double)(end-start)<<" Seconds"<<endl;
            Average_Syncrony << copling << '\t' << Total_syncrony_layer1 << '\t' << Total_syncrony_layer2 << endl;
            copling += data[6];//next copling
            Phases_layer2.close();
            Phases_layer1.close();
        }
        ofstream Last_Phase_layer1("Save/Last_Phase/layer1/Degree_Radian="+strDegree_Radian+".txt");
        ofstream Last_Phase_layer2("Save/Last_Phase/layer2/Degree_Radian="+strDegree_Radian+".txt");
        for (int i = 0; i < Number_of_node; i++) {
            Last_Phase_layer1 << Phases_layer1_previous[i] << endl;
            Last_Phase_layer2 << Phases_layer2_previous[i] << endl;
        }
        Average_Syncrony.close();
        Last_Phase_layer1.close();
        Last_Phase_layer2.close();
        delete Phases_layer1_previous;
        delete Phases_layer2_previous;
        delete Phases_next_layer1;
        delete Phases_next_layer2;
        delete Phases_history_delay_layer1[0];
        delete Phases_history_delay_layer2[0];
    }
    return 0;
}
