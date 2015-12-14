function [output] = analytic_sod(t)
%to solve Sod's Shock Tube problem
%reference: "http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html"
%   |       |   |     |         |
%   |       |   |     |         |
%   |       |   |     |         |
%___|_______|___|_____|_________|_______________
%   x1      x2  x0    x3        x4
%
%input require: t (time)
if nargin < 1
    %set default value
    t = 0.2;
end
%Initial conditions
x0 = 0.5;
rho_l = 1;
P_l = 1;
u_l = 0;

rho_r = 0.125;
P_r = 0.1;
u_r = 0;

gamma = 1.4;
g = sqrt( (gamma-1)/(gamma+1) );

%speeds of sound for left and right states
a_l = power( (gamma*P_l/rho_l),0.5);
a_r = power( (gamma*P_r/rho_r),0.5);

P_III = fzero('pressure_solve',pi);
v_cd = 2*(sqrt(gamma)/(gamma - 1))*(1 - power(P_III, (gamma - 1)/(2*gamma)));
rho_IV = rho_r*(( (P_III/P_r) + g^2 )/(1 + g*g*(P_III/P_r)));
v_shock = v_cd*((rho_IV/rho_r)/( (rho_IV/rho_r) - 1));
rho_III = (rho_l)*power((P_III/P_l),1/gamma);

%Key Positions
x1 = x0 - a_l*t;
x3 = x0 + v_cd*t;
x4 = x0 + v_shock*t;
%determining x2
c_2 = a_l - ((gamma - 1)/2)*v_cd;
x2 = x0 + (v_cd - c_2)*t;

%start setting values
n = 1000;   
%boundaries 
x_min = 0;
x_max = 1;

x = linspace(x_min,x_max,n);
output.x = x';
output.rho = zeros(n,1);   %density
output.P = zeros(n,1); %pressure
output.u = zeros(n,1); %velocity
output.e = zeros(n,1); %internal energy

for index = 1:n
    if output.x(index) < x1
        %Solution before x1, constant state
        output.rho(index) = rho_l;
        output.P(index) = P_l;
        output.u(index) = u_l;
    elseif (x1 <= output.x(index) && output.x(index) <= x2)
        %Solution in the zone between x1 and x2
        c = g*g*((x0 - output.x(index))/t) + (1 - g*g)*a_l; 
        output.rho(index) = rho_l*power((c/a_l),2/(gamma - 1));
        output.P(index) = P_l*power((output.rho(index)/rho_l),gamma);
        output.u(index) = (1 - g*g)*( (-(x0-output.x(index))/t) + a_l);
    elseif (x2 <= output.x(index) && output.x(index) <= x3)
        %Solution in the zone between x2 and x3
        output.rho(index) = rho_III;
        output.P(index) = P_III;
        output.u(index) = v_cd;
    elseif (x3 <= output.x(index) && output.x(index) <= x4)
        %Solution in the zone between x3 and x4
        output.rho(index) = rho_IV;
        output.P(index) = P_III;
        output.u(index) = v_cd;
    elseif x4 < output.x(index)
        %Solution after x4, constant state
        output.rho(index) = rho_r;
        output.P(index) = P_r;
        output.u(index) = u_r;
    end
    output.e(index) = output.P(index)/((gamma - 1)*output.rho(index));
end
end
