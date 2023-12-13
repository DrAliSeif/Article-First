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
| 0		| k_0=	| Coupling initial(frome exact initial to exact final was run)| 
| 0.01	| ∆k=		| Coupling step| 
| 3.0		| k_f=	| Coupling final| 
| 0.0		| τ_0= 	| delay start = number of data history| 
| 0.02	| ∆τ= 	| delay step| 
| 0.0		| τ_f= 	| delay end| 