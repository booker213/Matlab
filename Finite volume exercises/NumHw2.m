clear 
clc
format long

%% Define Constants

%Set Up grid
nc=100;  %Number of cells
nf=nc+1; %Number of edges
L=1;   %Length of domain
dx=L/nc; %Cell Length
x=0:dx:L; % Vector of edges
xc=dx/2:dx:L-dx/2;% Vector of centres

%Timestepping
dt=0.1;
