clc 
clear
ul=1.2; ur=0.4;
n=13; m=100
xi=linspace(-m,m,1000);
x=linspace(-m,m,1000);
u=zeros(1000,1);
figure
for t=0:.1:1
    for i=1:1000
if x(i)<-ul^n*t
    u(i)=ul;

elseif x(i)>-ur^n*t
u(i)=ur;

elseif -ul^n*t>x(i)>-ur^n*t
u(i)=(-x(i)/t)^(1/n);
end
    end
y= -ul^n*t
c=-ur^n*t
plot(x,u);
hold on
end
