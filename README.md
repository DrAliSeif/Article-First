# Article-First


### Calculate the information for two sine wave 

in file [Information/Calculate-sine-wave-information](https://github.com/DrAliSeif/Article-First/tree/main/Information/Calculate-sine-wave-information)

we have a two file .py
1. two-sine-wave-shifted.py (for create sine wave and plot it and save data in input_data)
2. calculate-information.py (for calculate each data in input_data file and put it in output_data) [in calculate-information.py code you can chose calculate all node or one of them for single
]

### create one layer network for kuramoto with cpp code

this code was parallel in c++
for run code you need ubuntu or server 
code for run: 	g++ main.cpp -fopenmp -o [name run]
				./[name run] &
				
				
## data.txt Syntax  				
example for data run

| Example        | Syntax      | Explane |
| ------|-----|-----|
| 1000	| N=		| Number of node| 
| 0		| t_0=	| First time| 
| 0.01	| ∆t=		| time step| 
| 20		| t_f=	| Final time| 
| 0		| k_0=	| Coupling initial (frome exact initial to exact final was run)| 
| 0.01	| ∆k=		| Coupling step| 
| 3.0		| k_f=	| Coupling final| 
| 0.0		| τ_0= 	| delay start = number of data history| 
| 0.02	| ∆τ= 	| delay step| 
| 0.0		| τ_f= 	| delay end| 

