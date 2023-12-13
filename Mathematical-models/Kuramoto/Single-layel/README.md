this code was parallel in c++
for run code you need ubuntu or server 
code for run: 	g++ main.cpp -fopenmp -o [name run]
				./[name run] &
				
				
## data.txt Syntax  				
example for data run

| Element        | Example        | Syntax      | Explane |
| ------|------|-----|-----|
| data[0]| 1000	| N=		| Number of node| 
| data[1]| 0		| t_0=	| First time| 
| data[2]| 0.01	| ∆t=		| time step| 
| data[3]| 20		| t_f=	| Final time| 
| data[4]| 0		| k_0=	| Coupling initial (frome exact initial to exact final was run)| 
| data[5]| 0.01	| ∆k=		| Coupling step| 
| data[6]| 3.0		| k_f=	| Coupling final| 
| data[7]| 0.0		| τ_0= 	| delay start = number of data history| 
| data[8]| 0.02	| ∆τ= 	| delay step| 
| data[9]| 0.0		| τ_f= 	| delay end| 