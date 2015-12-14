function y = pressure_solve(P)
%defines pressure function to be used in sod_tube
%Initial conditions defining left and right states
gamma = 1.4;

rho_l = 1;
P_l = 1;
u_l = 0;
rho_r = 0.125;
P_r = 0.1;
u_r = 0;



g = sqrt( (gamma-1)/(gamma+1) );

y = (P - P_r)*(( ((1 - g^2)^2)*((rho_r*(P + g*g*P_r))^-1) )^(0.5))...
    - 2*(sqrt(gamma)/(gamma - 1))*(1 - power(P, (gamma - 1)/(2*gamma)));
end
