#ifndef KURAMOTO_VERSION3_H_INCLUDED
#define KURAMOTO_VERSION3_H_INCLUDED
/************************************************************************************************/
/*** Topic: Kuramoto model with Runge-Kutta 4th Order Method Version 3                        ***/
/*** Date: 12/07/2022                                                                Ali-Seif ***/
/*** With Code:Blocks 20.03 and GCC compiler MinGW                                            ***/
/************************************************************************************************/
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#include<iostream>//for cout                                                                                              $$$$
#include<fstream>//infile /ofstream                                                                                       $$$$
#include <string>//for stod( )                                                                                            $$$$
#include <sstream>//stringstream ss(line)                                                                                 $$$$
#include<ctime>//For Example clock()                                                                                      $$$$
#include <cmath>//For Example pow                                                                                         $$$$
#include <omp.h>
#include <stdio.h>
//#define Pi 3.141592653589793238462643383279502884//pi number                                                              $$$$
using namespace std;//                                                                                                    $$$$
//------------------------------------------------------------------------------------------------------------------------$$$$
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                          Runge-Kutta 4th
//                                                          --------------
//                                                          \            /
//                                                           \          /
//                                                            \        /
//                                                             \      /
//                                                              \    /
//                                                               \  /
//                                                                \/
//-----------------------------------------------------------------------------------------------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
//@@@                                     dydt                                       @@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
double dydt(int Number_of_phase,                                                    //@@@                                   ---
            double frustration_inter_layer,                                         //@@@                                   ---
            int N,                                                                  //@@@                                   ---
            double dt,                                                              //@@@                                   ---
            double coupling,                                                        //@@@                                   ---
            double W,                                                               //@@@                                   ---
            const int* is_connected,                                                      //@@@                                   ---
            double phi,                                                             //@@@                                   ---
            double** phi_hist,                                                      //@@@                                   ---
            double secondlayer,
            double landa)                                                     //@@@                                   ---
{                                                                                   //@@@                                   ---
    double M = 0;                                                                   //@@@      connection calculated        ---
    double a = 0.0;                                                                 //@@@                                   ---
    for (int i = 0; i < N; i++){                                                    //@@@                                   ---
        a += (is_connected[i] * sin((phi_hist[i][0] - phi)));//               ---
    }                                                                               //@@@                                   ---
    // M = W + (coupling / (N * 1.0)) * a;                                          //@@@                                   ---
    //double landa = 10;                                                              //@@@                                   ---
    double connection = 0.0;                                                        //@@@                                   ---
    connection = (landa * sin(secondlayer - phi + frustration_inter_layer));        //@@@                                   ---
    M = W + (coupling / (N * 1.0)) * a + connection;                                //@@@               all sum             ---
    return M;                                                                       //@@@                                   ---
}                                                                                   //@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
//@@@                    Convert_next_to_history_and_previous                        @@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
void Convert_next_to_history_and_previous(int N,                                    //@@@                                   ---
                                          double dt,                                //@@@                                   ---
                                          double delay,                             //@@@                                   ---
                                          double** Phases_history_delay,            //@@@                                   ---
                                          double* Phases_next,                      //@@@                                   ---
                                          double* Phases_previous)                  //@@@                                   ---
{                                                                                   //@@@                                   ---
    int history = (delay / dt) + 1;                                                 //@@@                                   ---
    for (int i = 0; i < N; i++) {                                                   //@@@                                   ---
        for (int j = 0; j < history - 1; j++) {                                     //@@@                                   ---
            Phases_history_delay[i][j] = Phases_history_delay[i][j + 1];            //@@@                                   ---
        }                                                                           //@@@                                   ---
        Phases_history_delay[i][history - 1] = Phases_next[i];                      //@@@                                   ---
    }                                                                               //@@@                                   ---
    for (int i = 0; i < N; i++) { Phases_previous[i] = Phases_next[i]; }            //@@@                                   ---
}                                                                                   //@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
//@@@                                CCRK4                                           @@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
void Connected_Constant_Runge_Kutta_4(const double frustration_inter_layer,              //@@@                                   ---
                                      int N,                                        //@@@                                   ---
                                      double dt,                                    //@@@                                   ---
                                      double delay,                                 //@@@                                   ---
                                      double coupling,                              //@@@                                   ---
                                      const double* W,                                    //@@@                                   ---
                                      const int* const * adj,                                    //@@@                                   ---
                                      double* y,                                    //@@@                                   ---
                                      double** Phases_history_delay,                //@@@                                   ---
                                      double* Phases_secondlayer,                   //@@@                                   ---
                                      double* Phases_next,
                                      double landa)                          //@@@                                   ---
{                                                                                   //@@@                                   ---
    for (int i = 0; i < N; i++)                                                     //@@@                                   ---
    {                                                                               //@@@                                   ---
        double k1 = dydt(i, frustration_inter_layer, N,dt,//@                                   ---
                         coupling, W[i], adj[i], y[i], Phases_history_delay,        //@@@                                   ---
                         Phases_secondlayer[i],landa);                                    //@@@                                   ---
        double k2 = dydt(i, frustration_inter_layer, N,dt,//@                                   ---
                          coupling, W[i], adj[i], y[i] + k1 * dt / 2.0,             //@@@                                   ---
                          Phases_history_delay, Phases_secondlayer[i],landa);             //@@@                                   ---
        double k3 = dydt(i, frustration_inter_layer, N,dt,//@                                   ---
                          coupling, W[i], adj[i], y[i] + k2 * dt / 2.0,             //@@@                                   ---
                          Phases_history_delay, Phases_secondlayer[i],landa);             //@@@                                   ---
        double k4 = dydt(i, frustration_inter_layer, N,dt,//@                                   ---
                          coupling, W[i], adj[i], y[i] + k3 * dt, Phases_history_delay,//                                   ---
                           Phases_secondlayer[i],landa);                                  //@@@                                   ---
        Phases_next[i] = y[i] + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);         //@@@                                   ---
    }                                                                               //@@@                                   ---
    Convert_next_to_history_and_previous(N, dt, delay, Phases_history_delay,        //@@@                                   ---
                                         Phases_next, y);                           //@@@                                   ---
}                                                                                   //@@@                                   ---
//-----------------------------------------------------------------------------------------------------------------------------
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                            Scale_2_pi
//                                                          --------------
//                                                          \            /
//                                                           \          /
//                                                            \        /
//                                                             \      /
//                                                              \    /
//                                                               \  /
//                                                                \/
//-----------------------------------------------------------------------------------------------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
//@@@                                scale_2_pi phases                               @@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
void scale_2_pi(int N, double* phi)                                                 //@@@                                   ---
{                                                                                   //@@@                                   ---
    for (int i = 0; i < N; i++)                                                     //@@@                                   ---
    {                                                                               //@@@                                   ---
        phi[i] = phi[i] - 2 * M_PI * static_cast<int>(phi[i] / (2 * M_PI));             //@@@                                   ---
    }                                                                               //@@@                                   ---
}                                                                                   //@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
//@@@                                scale_pi phases                                 @@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
void scale_pi(int N, double* phi)                                                 //@@@                                   ---
{                                                                                   //@@@                                   ---
    for (int i = 0; i < N; i++)                                                     //@@@                                   ---
    {                                                                               //@@@                                   ---
        phi[i] -= 2 * M_PI * std::floor((phi[i] + M_PI) / (2 * M_PI));
    }                                                                               //@@@                                   ---

}                                                                                   //@@@                                   ---
//-----------------------------------------------------------------------------------------------------------------------------
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                         Order Parameter
//                                                          --------------
//                                                          \            /
//                                                           \          /
//                                                            \        /
//                                                             \      /
//                                                              \    /
//                                                               \  /
//                                                                \/
//-----------------------------------------------------------------------------------------------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
//@@@                                order_parameter                                 @@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
double order_parameter(int N, double* phi)                                          //@@@                                   ---
{                                                                                   //@@@                                   ---
    double rc = 0.0, rs = 0.0;                                                      //@@@                                   ---
    for (int j = 0; j < N; j++)                                                     //@@@                                   ---
    {                                                                               //@@@                                   ---
        rc += cos(phi[j]);                                                          //@@@                                   ---
        rs += sin(phi[j]);                                                          //@@@                                   ---
    }                                                                               //@@@                                   ---
    return sqrt(pow(rc, 2) + pow(rs, 2)) / (1.0 * N);                               //@@@                                   ---
}                                                                                   //@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
//@@@                             order_parameter_part                               @@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
double order_parameter_part(int Ini, int Fin, int N, double* phi)                   //@@@                                   ---
{                                                                                   //@@@                                   ---
    double rc = 0.0, rs = 0.0;                                                      //@@@                                   ---
    int conters = 0;                                                                //@@@                                   ---
    for (int j = (Ini - 1); j < (Fin); j++)                                         //@@@                                   ---
    {                                                                               //@@@                                   ---
        rc += cos(phi[j]);                                                          //@@@                                   ---
        rs += sin(phi[j]);                                                          //@@@                                   ---
        conters++;                                                                  //@@@                                   ---
    }                                                                               //@@@                                   ---
    return sqrt(pow(rc, 2) + pow(rs, 2)) / (1.0 * conters);                         //@@@                                   ---
}                                                                                   //@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
//@@@                             order_parameter_specific                           @@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
double order_parameter_specific(double Ini, double Fin, int N, double* phi, double* frequency)//                            ---
{                                                                                   //@@@                                   ---
    double rc = 0.0, rs = 0.0;                                                      //@@@                                   ---
    int conters = 0;                                                                //@@@                                   ---
    for (int j = 0; j < N; j++)                                                     //@@@                                   ---
    {                                                                               //@@@                                   ---
        if (frequency[j] <= Ini && Fin <= frequency[j]) {                           //@@@                                   ---
            rc += cos(phi[j]);                                                      //@@@                                   ---
            rs += sin(phi[j]);                                                      //@@@                                   ---
            conters++;                                                              //@@@                                   ---
        }                                                                           //@@@                                   ---
    }                                                                               //@@@                                   ---
    return sqrt(pow(rc, 2) + pow(rs, 2)) / (1.0 * conters);                         //@@@                                   ---
}                                                                                   //@@@                                   ---
//-----------------------------------------------------------------------------------------------------------------------------
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                   Hystory and previous phases
//                                                          --------------
//                                                          \            /
//                                                           \          /
//                                                            \        /
//                                                             \      /
//                                                              \    /
//                                                               \  /
//                                                                \/
//-----------------------------------------------------------------------------------------------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@Create history for delay           ---
//@@@                             Create Hystory phases                              @@@@Create history for delay           ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@(Phases history delay & previous)  ---
double** Hystorydelay_phases(const int N, const int history, const double* Phases_initial) {          //@@@calculate initial theta            ---
    double** Phases_history_delay = new double* [N];                                //@@@                                   ---
    for (int i = 0; i < N; i++) {                                                   //@@@                                   ---
        Phases_history_delay[i] = new double[history];                              //@@@[node][delay]                      ---
    }                                                                               //@@@                                   ---
    for (int i = 0; i < N; i++) {                                                   //@@@                                   ---
        Phases_history_delay[i][0] = Phases_initial[i];                             //@@@                                   ---
    }                                                                               //@@@                                   ---
    for (int t = 1; t < history; t++)                                               //@@@                                   ---
    {                                                                               //@@@                                   ---
        for (int i = 0; i < N; i++)                                                 //@@@                                   ---
        {                                                                           //@@@                                   ---
            if (Phases_history_delay[i][t - 1] >= (double)(M_PI / 2.0))               //@@@                                   ---
            {                                                                       //@@@                                   ---
                Phases_history_delay[i][t] = Phases_history_delay[i][t - 1] - (double)(1.5 * M_PI);//                         ---
            }                                                                       //@@@                                   ---
            else                                                                    //@@@                                   ---
            {                                                                       //@@@                                   ---
                Phases_history_delay[i][t] = Phases_history_delay[i][t - 1] + (double)(M_PI / 2.0);//                         ---
            }                                                                       //@@@                                   ---
        }                                                                           //@@@                                   ---
    }                                                                               //@@@                                   ---
    return Phases_history_delay;                                                    //@@@                                   ---
}                                                                                   //@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
//@@@                                  previous phases                               @@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
double* previous_phases(int N, int history, double** Phases_history_delay) {        //@@@calculate initial theta            ---
    double* previous = new double[N];                                               //@@@Making Array                       ---
    for (int i = 0; i < N; i++)                                                     //@@@                                   ---
    {                                                                               //@@@                                   ---
        if (Phases_history_delay[i][history - 1] >= (double)(M_PI / 2.0))             //@@@                                   ---
        {                                                                           //@@@                                   ---
            previous[i] = Phases_history_delay[i][history - 1] - (double)(1.5 * M_PI);//@@@                                   ---
        }                                                                           //@@@                                   ---
        else                                                                        //@@@                                   ---
        {                                                                           //@@@                                   ---
            previous[i] = Phases_history_delay[i][history - 1] + (double)(M_PI / 2.0);//@@@                                   ---
        }                                                                           //@@@                                   ---
    }                                                                               //@@@                                   ---
                                                                                    //@@@                                   ---
    return previous;                                                                //@@@                                   ---
}                                                                                   //@@@                                   ---
//-----------------------------------------------------------------------------------------------------------------------------
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                        read file to arrey
//                                                          --------------
//                                                          \            /
//                                                           \          /
//                                                            \        /
//                                                             \      /
//                                                              \    /
//                                                               \  /
//                                                                \/
//-----------------------------------------------------------------------------------------------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
//@@@                               count rows file in .txt                          @@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Read data from text               ---
double* read_initial_1D(string Filename, int Numberofnode)                          //@@@ (Phases & frequency & Matrix)     ---
{                                                                                   //@@@                                   ---
    double* data_1D = new double[Numberofnode];                                     //@@@                                   ---
    ifstream file("Example/" + Filename + ".txt");                                  //@@@                                   ---
    for (int i = 0; i < Numberofnode; i++)                                          //@@@                                   ---
    {                                                                               //@@@                                   ---
        file >> data_1D[i];                                                         //@@@                                   ---
    }                                                                               //@@@                                   ---
    //cout << Filename + "\t\tloaded \t core number= "<<omp_get_thread_num()<< endl;                                        //@@@                                   ---
    file.close();                                                                   //@@@                                   ---
    return data_1D;                                                                 //@@@                                   ---
}                                                                                   //@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
//@@@                          Read matrix connection                               //@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
int** read_initial_2D(string Filename, int Numberofnode)                            //@@@                                   ---
{                                                                                   //@@@                                   ---
    int** data_2D = new int* [Numberofnode];                                        //@@@                                   ---
    for (int i = 0; i < Numberofnode; i++)                                          //@@@                                   ---
        data_2D[i] = new int[Numberofnode];                                         //@@@                                   ---
    ifstream file("Example/" + Filename + ".txt");                                  //@@@                                   ---
    if (!file)                                                                      //@@@                                   ---
    {                                                                               //@@@                                   ---
        cout << "Matrix file is not here!" << endl;                                 //@@@                                   ---
        return data_2D;                                                             //@@@                                   ---
    }                                                                               //@@@                                   ---
    else                                                                            //@@@                                   ---
    {                                                                               //@@@                                   ---
        for (int i = 0; i < Numberofnode; i++)                                      //@@@                                   ---
        {                                                                           //@@@                                   ---
            for (int j = 0; j < Numberofnode; j++)                                  //@@@                                   ---
            {                                                                       //@@@                                   ---
                int elem = 0;                                                       //@@@                                   ---
                file >> elem;                                                       //@@@                                   ---
                data_2D[i][j] = elem;                                               //@@@                                   ---
            }                                                                       //@@@                                   ---
        }                                                                           //@@@                                   ---
    }                                                                               //@@@                                   ---
    cout << Filename + "\t\tloaded" << endl;                                        //@@@                                   ---
    return data_2D;                                                                 //@@@                                   ---
}                                                                                   //@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
//@@@                          Read matrix frustration                              //@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
double** read_initial_2D_double(string Filename, int Numberofnode)                  //@@@                                   ---
{                                                                                   //@@@                                   ---
    double** data_2D = new double* [Numberofnode];                                  //@@@                                   ---
    for (int i = 0; i < Numberofnode; i++)                                          //@@@                                   ---
        data_2D[i] = new double[Numberofnode];                                      //@@@                                   ---
    ifstream file("Example/" + Filename + ".txt");                                  //@@@                                   ---
    if (!file)                                                                      //@@@                                   ---
    {                                                                               //@@@                                   ---
        cout << "Matrix file is not here!" << endl;                                 //@@@                                   ---
        return data_2D;                                                             //@@@                                   ---
    }                                                                               //@@@                                   ---
    else                                                                            //@@@                                   ---
    {                                                                               //@@@                                   ---
        for (int i = 0; i < Numberofnode; i++)                                      //@@@                                   ---
        {                                                                           //@@@                                   ---
            for (int j = 0; j < Numberofnode; j++)                                  //@@@                                   ---
            {                                                                       //@@@                                   ---
                double elem = 0.0;                                                  //@@@                                   ---
                file >> elem;                                                       //@@@                                   ---
                data_2D[i][j] = elem;                                               //@@@                                   ---
            }                                                                       //@@@                                   ---
        }                                                                           //@@@                                   ---
    }                                                                               //@@@                                   ---
    cout << Filename + "\t\tloaded" << endl;                                        //@@@                                   ---
    return data_2D;                                                                 //@@@                                   ---
}                                                                                   //@@@                                   ---
//-----------------------------------------------------------------------------------------------------------------------------
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                             data.txt
//                                                          --------------
//                                                          \            /
//                                                           \          /
//                                                            \        /
//                                                             \      /
//                                                              \    /
//                                                               \  /
//                                                                \/
//
//-----------------------------------------------------------------------------------------------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
//@@@                               count rows file in .txt                         //@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
int count_rows_file(string file1)                                                   //@@@                                   ---
{                                                                                   //@@@                                   ---
    int rows = 0, cols = 0;                                                         //@@@                                   ---
    string line, item;                                                              //@@@                                   ---
    ifstream file(file1);                                                           //@@@                                   ---
    while (getline(file, line))                                                     //@@@                                   ---
    {                                                                               //@@@                                   ---
        rows++;                                                                     //@@@                                   ---
        if (rows == 1)                                                              //@@@First row only:                    ---
        {                                                                           //@@@determine the number of columns    ---
            stringstream ss(line);                                                  //@@@Set up up a stream from this line  ---
            while (ss >> item) cols++;                                              //@@@Each item delineated by spaces     ---
        }                                                                           //@@@adds one to cols                   ---
    }                                                                               //@@@                                   ---
    file.close();                                                                   //@@@                                   ---
    cout << "\nFile had " << rows << " rows and " << cols << " columns\n" << endl;  //@@@                                   ---
    return rows;                                                                    //@@@                                   ---
}                                                                                   //@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
//@@@                               Read Data in .txt                               //@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
double* read_data(double* data)                                                     //@@@                                   ---
{                                                                                   //@@@                                   ---
    string kk;                                                                      //@@@                                   ---
    ifstream fp("data.txt");                                                        //@@@                                   ---
    if (!fp)                                                                        //@@@                                   ---
    {                                                                               //@@@                                   ---
        cout << "Data file is not here!" << endl;                                   //@@@                                   ---
    }                                                                               //@@@                                   ---
    else                                                                            //@@@                                   ---
    {                                                                               //@@@                                   ---
        fp >> kk;                                                                   //@@@                                   ---
        string line, item;                                                          //@@@                                   ---
        int i = 1;                                                                  //@@@                                   ---
        while (getline(fp, line))                                                   //@@@                                   ---
        {                                                                           //@@@                                   ---
            fp >> kk;                                                               //@@@                                   ---
            data[i] = stod(kk);                                                     //@@@                                   ---
            i++;                                                                    //@@@                                   ---
        }                                                                           //@@@                                   ---
    }                                                                               //@@@                                   ---
    cout << "|----------------------------------------------------------------------------------------------------|" << endl;
    cout << "|Data=>>                              \t\t\t\t\t\t\t\t     |" << endl;      //                                 ---
    cout << "|                                      \t\t\t\t\t\t\t\t     |" << endl;     //                                 ---
    cout << "|data[1]= Number =\t\t" << data[1] <<"\t\t\t\t\t\t\t\t     |" << endl;      //                                 ---
    cout << "|data[2]= First Time =\t\t" << data[2] <<"\t\t\t\t\t\t\t\t     |" << endl;  //                                 ---
    cout << "|data[3]= dt =\t\t\t" << data[3] <<"\t\t\t\t\t\t\t\t     |" << endl;        //                                 ---
    cout << "|data[4]= Final Time =\t\t" << data[4] << "\t\t\t\t\t\t\t\t     |" << endl; //                                 ---
    cout << "|----------------------------------------------------------------------------------------------------|" << endl;
    cout << "|coupling=>>                              \t\t\t\t\t\t\t     |" << endl;    //                                 ---
    cout << "|                                      \t\t\t\t\t\t\t\t     |" << endl;     //                                 ---
    cout << "|data[5]= coupling_start =\t" << data[5] <<"\t\t\t\t\t\t\t\t     |" << endl;//                                 ---
    cout << "|data[6]= coupling_step =\t" << data[6] <<"\t\t\t\t\t\t\t\t     |" << endl; //                                 ---
    cout << "|data[7]= coupling_end =\t" << data[7] <<"\t\t\t\t\t\t\t\t     |" << endl;  //                                 ---
    cout << "|----------------------------------------------------------------------------------------------------|" << endl;
    cout << "|delay=>>                              \t\t\t\t\t\t\t\t     |" << endl;     //                                 ---
    cout << "|                                      \t\t\t\t\t\t\t\t     |" << endl;     //                                 ---
    cout << "|data[8]= delay_start =\t\t" << data[8] << "\t\t\t\t\t\t\t\t     |" << endl;//                                 ---
    cout << "|data[9]= delay_step =\t\t" << data[9] << "\t\t\t\t\t\t\t\t     |" << endl; //                                 ---
    cout << "|data[10]= delay_end =\t\t" << data[10] << "\t\t\t\t\t\t\t\t     |" << endl;//                                 ---
    cout << "|----------------------------------------------------------------------------------------------------|" << endl;
    fp.close();                                                                     //@@@                                   ---
    return data;                                                                    //@@@                                   ---
}                                                                                   //@@@                                   ---
//-----------------------------------------------------------------------------------------------------------------------------
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                              |    |
//                                                            Correlation
//                                                          --------------
//                                                          \            /
//                                                           \          /
//                                                            \        /
//                                                             \      /
//                                                              \    /
//                                                               \  /
//                                                                \/
//-----------------------------------------------------------------------------------------------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
//@@@                              calcoulate correlation                           //@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
double** Correlation(int N, double* Phases) {                                       //@@@                                   ---
    double** Correlation_Matrix = new double* [N];                                  //@@@                                   ---
    for (int i = 0; i < N; i++)                                                     //@@@                                   ---
        Correlation_Matrix[i] = new double[N];                                      //@@@                                   ---
    for (int i = 0; i < N; i++) {                                                   //@@@                                   ---
        for (int j = 0; j < N; j++) {                                               //@@@                                   ---
            Correlation_Matrix[i][j] = abs(cos(Phases[i] - Phases[j]));             //@@@                                   ---
        }                                                                           //@@@                                   ---
    }                                                                               //@@@                                   ---
    return Correlation_Matrix;                                                      //@@@                                   ---
    delete[] Correlation_Matrix;                                                    //@@@                                   ---
}                                                                                   //@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
//@@@                               print it                                        //@@@                                   ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                   ---
void Print_2D(string stringname, int N, string strcopling, string strdelay,         //@@@                                   ---
               double** data) {                                                     //@@@                                   ---
    ofstream Print("Save/Correlation/delay=" + strdelay + "_copling=" + strcopling +//@@@                                   ---
                    stringname +"(Node=)VS(Node=).txt");                            //@@@                                   ---
    for (int i = 0; i < N; i++) {                                                   //@@@                                   ---
        for (int j = 0; j < N; j++) {                                               //@@@                                   ---
            Print << data[i][j] << '\t';                                            //@@@                                   ---
        }                                                                           //@@@                                   ---
        Print << endl;                                                              //@@@                                   ---
    }                                                                               //@@@                                   ---
}                                                                                   //@@@                                   ---
//-----------------------------------------------------------------------------------------------------------------------------
#endif // KURAMOTO_VERSION3_H_INCLUDED
