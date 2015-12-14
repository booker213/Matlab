clear all

%%Define IC and constants
J=20; dx=1/J; dt=0.001;
s=(dt^2)/(dx^2);
c=4;
steps=1125;
time=steps*dt;
k=c*s;

x=linspace(0+dx,1-dx,J-1);
u0=zeros(J-1,1);
for i=1:J-1
u0(i)=sin(2*pi*x(i));
end


Uprev=u0;
U=zeros(J-1,1);
%%scheme to initialise u(i=1,j)
U(1)=u0(1)+(k/2)*(u0(2)-2*u0(1));          %%u0(i-1)=0
U(J-1)=u0(J-1)+(k/2)*(-2*u0(J-1)+u0(J-2));  %%u0(J)=0
for i=2:J-2
U(i)=u0(i)+(k/2)*(u0(i+1)-2*u0(i)+u0(i-1));
end

A=zeros(J-1,J-1);
%%row 0 and row j+1 are zeroes bc BCs define u0=uj=0 for all t
A(1,1)=2-2*k;
A(1,2)=k;
A(J-1,J-2)=k;
A(J-1,J-1)=2-2*k;

for j=2:J-2
A(j,j-1)=k;
A(j,j+1)=k;
A(j,j)=2-2*k;
end

Ustorage=zeros(J-1,1);
for n=1:steps-1
Ustorage=U;
U=A*U-Uprev;
Uprev=Ustorage;
end


plot(x,U,'b');
title('time=1.25')
xlabel('x')
ylabel('u{x,t)')
hold on

