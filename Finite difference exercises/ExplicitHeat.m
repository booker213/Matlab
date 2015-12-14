clear all

%%Define IC and constants
J=20; dx=1/J; dt=0.0012;
mu=dt/(dx^2);
steps=10;
time=steps*dt;

x=linspace(0,1,J+1);
u0=zeros(J+1,1);
for i=1:21
if x(i)<0.5
   u0(i)=2*x(i);
elseif x(i)>=0.5
   u0(i)=2-2*x(i);
end
end

U=u0;
A=zeros(J+1,J+1);
%%row 0 and row j+1 are zeroes bc BCs define u0=uj=0 for all t
A(2,2)=1-2*mu;
A(2,3)=mu;
A(J,J-1)=mu;
A(J,J)=1-2*mu;

for j=3:J-1
A(j,j-1)=mu;
A(j,j+1)=mu;
A(j,j)=1-2*mu;
end

for n=1:steps
U=A*U;
end


plot(x,U,'b');
hold on
uexact=zeros(1,J+1);
for i=1:21
syms n;
uexact(i)= symsum(exp(-((n*pi)^2)*time)*((8*sin(n*pi/2))/(pi*n)^2)*sin(n*pi*x(i)),1,50);
end

plot(x,uexact,'r');
