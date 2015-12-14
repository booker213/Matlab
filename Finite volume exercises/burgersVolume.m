%finite volume for burgers equation, only one speed...U(x,t)
%doing this to then adapt for homework. 

%order - parameter values, define mesh, compute cell areas and side
%vectors, inital values on mesh, insert ghost cells, (*) start time loop,
%ghost cell values, calculate time step, calculate inner fluxes, implement
%solver and update solutions, return to (*) til end of time, output results

%parameters and defining mesh

tend = 2; %end of time, two seconds
L=2*pi; %domain length [0,1]
nt = 1000; %number of time steps, t_n = 1,...10
xk = 100; %number of cells (center of cells)


cells = 1:xk; %vector numbering each cell. may not be needed
dx = L/xk; %for simple mesh where diff in x is the same
celledge = 0:dx:L; %vector of values at the edge of the cells
cellcenters = zeros(1,xk); %vectors of cell centre value
for i = 1:xk
    cellcenters(i) = (celledge(i+1)+celledge(i))/2;
end    

dt = tend/nt;

cellareas = zeros(1,xk); %vector for area of diff cells. Long winded as all same in this simple mesh
for i =1:xk
    cellareas(i) = (celledge(i+1)-celledge(i))*dt;
end


%initial condition sin(x)

u0 = zeros(1,xk);
for i=1:xk
    u0(i)= sin(cellcenters(i));
end
Unew=u0;
U=u0;

for q = 1:nt
    
    for i = 2:xk-2
        
        Unew(i) = U(i) - (dt/dx)*(reimann(U(i),U(i+1))-reimann(U(i-1),U(i)));
        
    end
    U = Unew;
    plot(cells,U) 
    
    
end


