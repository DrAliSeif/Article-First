// C++ Program to demonstrate Mathematical model (kuramoto double layer)
#include"Kuramoto.Version4.h"//library Kuramoto version 4 (ubuntu version push in github)

int main() {
    // Hint1: count_rows_cols_file: para in address of file that is ./data.txt
    // Hint2: read_data: first para is number of rows in data file and second para is boolean[0=dont show data,1=show data]
    // data[0]=N & data[1]=L & data[2]=a
    // data[3]=t_0 & data[4]=∆t & data[5]=t_f
    // data[6]=k_0 & data[7]=∆k & data[8]=k_f
    // data[9]=τ_0 & data[10]=∆τ & data[11]=τ_f
    double* data=read_data(count_rows_cols_file("data.txt"),0);
    double* frequency_layer1 = read_initial_1D("W=Natural frequency/Layer1_mean1", int(data[0]));
    double* frequency_layer2 = read_initial_1D("W=Natural frequency/Sarika0.8a_mean1", int(data[0]));

    double* Phases_initial_layer1 = read_initial_1D("P=Initial Phases/Phases_initial_layer1_origin", int(data[0]));//Initial Phases  P
    double* Phases_initial_layer2 = read_initial_1D("P=Initial Phases/Phases_initial_layer2_origin", int(data[0]));//Initial Phases  P
    int** adj_layer1 = read_initial_2D("A=Intralayer adjacency matrix/Layer1", int(data[0]));//adjacency matrix  A
    int** adj_layer2 = read_initial_2D("A=Intralayer adjacency matrix/Layer2", int(data[0]));//adjacency matrix  A
    double Delay_variable = data[9];
    // while (Delay_variable < (data[11])){Delay_variable+=data[10]} // Delay loop
    double* Phases_next_layer1 = new double[int(data[0])];
    double* Phases_next_layer2 = new double[int(data[0])];
    double** Phases_history_delay_layer1 = memory_of_delay_of_phases(int(data[0]),Delay_variable,data[4],Phases_initial_layer1);
    double** Phases_history_delay_layer2 = memory_of_delay_of_phases(int(data[0]),Delay_variable,data[4],Phases_initial_layer2);
    double* Phases_layer1_previous = shift_pi2_phases(int(data[0]),Delay_variable,data[4], Phases_history_delay_layer1);//Phases changer
    double* Phases_layer2_previous = shift_pi2_phases(int(data[0]),Delay_variable,data[4], Phases_history_delay_layer2);//Phases changer
    // Hint3: When i change it that add variable to data.txt 
    ofstream Avg_Sync(name_file_data("Save/Avg_Sync/",data,12)+".txt");
    double Coupling_variable = data[6];
    cout<<"8. G to Coupling_variable. :)"<<endl;
    while (Coupling_variable <= (data[8])) { // Coupling loop
        ofstream Save_phases1_for_each_coupling("Save/Phases/layer1/k="+to_string(Coupling_variable)+".txt");
        ofstream Save_phases2_for_each_coupling("Save/Phases/layer2/k="+to_string(Coupling_variable)+".txt");
        double Total_synchrony_layer1 = 0;
        double Total_synchrony_layer2 = 0;
        int counter_of_total_sync =0;
        double Time_variable = data[3];// reset time for new time
        while (Time_variable < (data[5] + data[4])) {
            Connected_Constant_Runge_Kutta_4(data,Delay_variable, Coupling_variable, frequency_layer1, adj_layer1, Phases_layer1_previous, Phases_history_delay_layer1,Phases_next_layer1, Phases_layer2_previous);
            Connected_Constant_Runge_Kutta_4(data,Delay_variable, Coupling_variable, frequency_layer2, adj_layer2, Phases_layer2_previous, Phases_history_delay_layer2,Phases_next_layer2, Phases_layer1_previous);
            double synchrony_layer1 = order_parameter(int(data[0]), Phases_layer1_previous);// order parameters
            double synchrony_layer2 = order_parameter(int(data[0]), Phases_layer2_previous);// order parameters
            Save_phases1_for_each_coupling << Time_variable << '\t';
            Save_phases2_for_each_coupling << Time_variable << '\t';
            for (int i = 0; i < int(data[0]); i++) {
                Save_phases1_for_each_coupling << Phases_layer1_previous[i] << '\t';
                Save_phases2_for_each_coupling << Phases_layer2_previous[i] << '\t';
            }
            Save_phases1_for_each_coupling << endl;
            Save_phases2_for_each_coupling << endl;
            if (Time_variable >= int(data[5] * 0.8)) {// add sync to total sync
                Total_synchrony_layer1 += synchrony_layer1;
                Total_synchrony_layer2 += synchrony_layer2;
                counter_of_total_sync+=1;
            }
            Time_variable += data[4];
        }
        Total_synchrony_layer1=Total_synchrony_layer1/counter_of_total_sync;// calculate total sync and pint it
        Total_synchrony_layer2=Total_synchrony_layer2/counter_of_total_sync;// calculate total sync and pint it
        cout<< Coupling_variable << '\t' << Total_synchrony_layer1 << '\t' << Total_synchrony_layer2<<endl;
        Avg_Sync << Coupling_variable << '\t' << Total_synchrony_layer1<< '\t' << Total_synchrony_layer2<< endl;
        Save_phases1_for_each_coupling.close();
        Save_phases2_for_each_coupling.close();
        Coupling_variable += data[7];// next Coupling_variable
    }
    write_last_phase("Save/Last_Phase/layer1/",data,12,Phases_layer1_previous);
    write_last_phase("Save/Last_Phase/layer2/",data,12,Phases_layer2_previous);
    Avg_Sync.close();
    delete Phases_layer1_previous;
    delete Phases_next_layer1;
    delete Phases_history_delay_layer1[0];
    delete Phases_layer2_previous;
    delete Phases_next_layer2;
    delete Phases_history_delay_layer2[0];
    return 0;
}
