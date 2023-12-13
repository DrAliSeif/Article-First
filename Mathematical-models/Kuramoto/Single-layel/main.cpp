// C++ Program to demonstrate Mathematical model (kuramoto single layer)
#include"Kuramoto.Version4.h"//library Kuramoto version 4 (ubuntu version push in github)
//#include <time.h>

#include<iostream>//for cout
using namespace std;


int main() {
    // Read data from data.txt and write them to pointer 1d
    // Hint1: count_rows_cols_file: para in address of file that is ./data.txt
    // Hint2: read_data: first para is number of rows in data file and second para is boolean[0=dont show data,1=show data]
    // data[0]=N & data[1]=L & data[2]=a
    // data[3]=t_0 & data[4]=∆t & data[5]=t_f
    // data[6]=k_0 & data[7]=∆k & data[8]=k_f
    // data[9]=τ_0 & data[10]=∆τ & data[11]=τ_f
    double* data=read_data(count_rows_cols_file("data.txt"),0);


    double* frequency_layer1 = read_initial_1D("W=Natural frequency/Layer1", int(data[0]));
    double* Phases_initial_layer1 = read_initial_1D("P=Initial Phases/Phases_initial_layer1_origin", int(data[0]));//Initial Phases  P
    int** adj_layer1 = read_initial_2D("A=Intralayer adjacency matrix/Layer1", int(data[0]));//adjacency matrix  A
    double Delay_variable = data[9];
    // while (Delay_variable < (data[11])){Delay_variable+=data[10]} // Delay loop
    double* Phases_next_layer1 = new double[int(data[0])];
    double** Phases_history_delay_layer1 = memory_of_delay_of_phases(int(data[0]),Delay_variable,data[4],Phases_initial_layer1);

    //double* Phases_layer1_previous = previous_phases(int(data[0]), (int(Delay_variable / data[4]) + 1), Phases_history_delay_layer1);//Phases changer
   
    ofstream Avg_Sync(name_file_data(data,12));
    double Coupling_variable = data[6];
    /*   while (Coupling_variable < (data[8])) { // Coupling loop
        time_t start = time(NULL);
        ostringstream ostrcopling;// declaring output string stream
        ostrcopling << Coupling_variable;// Sending a number as a stream output
        string strcopling = ostrcopling.str();// the str() converts number into string
        ofstream Phases_layer2("Save/Phases/Degree_Radian="+strDegree_Radian+"_copling="+strcopling+"layer2(time)VS(Node).txt");//     ---
        ofstream Phases_layer1("Save/Phases/Degree_Radian="+strDegree_Radian+"_copling="+strcopling+"layer1(time)VS(Node).txt");
        double Total_syncrony_layer1 = 0;
        double Total_syncrony_layer2 = 0;
        // time loop
        double time_loop = double(data[3]);// reset time for new time
        while (time_loop < (data[5] + data[4])) {
            // Runge Kutta and scaling phase
            Connected_Constant_Runge_Kutta_4(data[2], int(data[0]), data[4], Delay_variable, Coupling_variable, frequency_layer1, adj_layer1, Phases_layer1_previous, Phases_history_delay_layer1, Phases_layer2_previous,Phases_next_layer1,data[1]);
            Connected_Constant_Runge_Kutta_4(data[2], int(data[0]), data[4], Delay_variable, Coupling_variable, frequency_layer2, adj_layer2, Phases_layer2_previous, Phases_history_delay_layer2, Phases_layer1_previous,Phases_next_layer2,data[1]);
            scale_pi(int(data[0]), Phases_layer1_previous);// change values to be in range 0 to 2*Pi
            scale_pi(int(data[0]), Phases_layer2_previous);
            double syncrony_layer1 = order_parameter(int(data[0]), Phases_layer1_previous);// order parameters
            double syncrony_layer2 = order_parameter(int(data[0]), Phases_layer2_previous);// order parameters
            Phases_layer2 << time_loop << '\t';//Phases_layer1 << time << '\t';
            Phases_layer1 << time_loop << '\t';
            for (int i = 0; i < int(data[0]); i++) {
                Phases_layer2 << Phases_layer2_previous[i] << '\t';//Phases_layer1 << Phases_layer1_previous[i] << '\t';
                Phases_layer1 << Phases_layer1_previous[i] << '\t';
            }
            Phases_layer2 << endl;// Phases_layer1 << endl;
            Phases_layer1 << endl;
            if (time_loop > int(data[5] * 0.2)) {// add sync to total sync
                Total_syncrony_layer1 += syncrony_layer1;
                Total_syncrony_layer2 += syncrony_layer2;
            }
            time_loop += data[4];
        }
        Total_syncrony_layer1=Total_syncrony_layer1/(int(int(data[5] * 0.2) / data[4])*4);// calculate total sync and pint it
        Total_syncrony_layer2=Total_syncrony_layer2/(int(int(data[5] * 0.2) / data[4])*4);
        time_t end = time(NULL);
        cout<< Coupling_variable << '\t' << Total_syncrony_layer2 <<'\t' <<"Execution Time: "<< (double)(end-start)<<" Seconds"<<endl;
        Avg_Sync << Coupling_variable << '\t' << Total_syncrony_layer1 << '\t' << Total_syncrony_layer2 << endl;
        Coupling_variable += data[7];//next Coupling_variable
        Phases_layer2.close();
        Phases_layer1.close();
    }
    ofstream Last_Phase_layer1("Save/Last_Phase/layer1/Degree_Radian="+strDegree_Radian+".txt");
    ofstream Last_Phase_layer2("Save/Last_Phase/layer2/Degree_Radian="+strDegree_Radian+".txt");
    for (int i = 0; i < int(data[0]); i++) {
        Last_Phase_layer1 << Phases_layer1_previous[i] << endl;
        Last_Phase_layer2 << Phases_layer2_previous[i] << endl;
    }
    Avg_Sync.close();
    Last_Phase_layer1.close();
    Last_Phase_layer2.close();
    delete Phases_layer1_previous;
    delete Phases_layer2_previous;
    delete Phases_next_layer1;
    delete Phases_next_layer2;
    delete Phases_history_delay_layer1[0];
    delete Phases_history_delay_layer2[0];*/
    return 0;
}
