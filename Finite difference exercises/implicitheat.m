clear all

%%Define IC and constants
J=20; dx=1/J; dt=0.001;
mu=dt/(dx^2);
time=0.5;
steps=time/dt;


x=linspace(0+dx,1-dx,J-1);
u0=zeros(J-1,1);
for i=1:J-1
u0(i)=100*sin(pi*x(i));
end

U=u0;
A=zeros(J-1,J-1);
%%row 0 and row j+1 are zeroes bc BCs define u0=uj=0 for all t
A(1,1)=1+2*mu;
A(1,2)=-mu;
A(J-1,J-2)=-mu;
A(J-1,J-1)=1+2*mu;

for j=2:J-2
A(j,j-1)=-mu;
A(j,j+1)=-mu;
A(j,j)=1+2*mu;
end

for n=1:steps
U=A\U;
end


plot(x,U,'b');
hold on

uexact=zeros(1,J-1);
for i=1:J-1
uexact(i)= 100*exp(-((pi)^2)*time)*sin(pi*x(i));
end

plot(x,uexact,'r');