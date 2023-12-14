# Article-First


### Calculate the information for two sine wave 

in file [Information/Calculate-sine-wave-information](https://github.com/DrAliSeif/Article-First/tree/main/Information/Calculate-sine-wave-information)

we have a two file .py
1. two-sine-wave-shifted.py (for create sine wave and plot it and save data in input_data)
![](https://github.com/DrAliSeif/Article-First/blob/main/Information/Calculate-sine-wave-information/input_data/two_sine_shifted/two_sine_wave_with_shifted_pi2_noise.png).
2. calculate-information.py (for calculate each data in input_data file and put it in output_data) [in calculate-information.py code you can chose calculate all node or one of them for single
]
![](https://github.com/DrAliSeif/Article-First/blob/main/Information/Calculate-sine-wave-information/output_data/10ColsRandomGaussian-1.png).


### create one layer network for kuramoto with cpp code

this code was parallel in c++
for run code you need ubuntu or server 
code for run: 	g++ main.cpp -fopenmp -o [name run]
				./[name run] &
				
				
## data.txt Syntax  	


1. Initial phases is in [-Pi,Pi] in file [input_data/P=Initial Phases/Phases_initial_layer1_origin.txt](https://github.com/DrAliSeif/Article-First/tree/main/Mathematical-models/Kuramoto/Single-layel/input_data/P%3DInitial%20Phases)


example for data run

| Element        | Example        | Syntax      | Explain |
| ------|------|-----|-----|
| data[0]| 1000	| N=		| Number of node| 
| data[1]| 10	| L=		| Coupling inter layer| 
| data[2]| 1.57	| a=		| Frustration| 
| data[3]| 0		| t_0=	| First time| 
| data[4]| 0.01	| ∆t=		| time step| 
| data[5]| 20		| t_f=	| Final time| 
| data[6]| 0		| k_0=	| Coupling initial (from exact initial to exact final was run)| 
| data[7]| 0.01	| ∆k=		| Coupling step| 
| data[8]| 3.0		| k_f=	| Coupling final| 
| data[9]| 0.0		| τ_0= 	| delay start = number of data history| 
| data[10]| 0.02	| ∆τ= 	| delay step| 
| data[11]| 0.0		| τ_f= 	| delay end| 
