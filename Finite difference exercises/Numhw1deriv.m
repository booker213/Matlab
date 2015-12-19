function y = Numhw1deriv(B)
%% Constants
q=0.99; alpha=0.4709; beta = 1;
H=1; Bbase= 1.178164343; Bt=0.585373798;


%% Derivative

y= (alpha*B^3-q)/(beta*B^3);
end