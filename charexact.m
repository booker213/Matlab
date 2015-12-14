clc , clear;

xi=linspace(0,pi,2000);

x=zeros(2000,1);
u=zeros(2000,1);
u2=zeros(2000,1);
u3=zeros(2000,1);
for i=1:2000
if xi(i)<1
x(i)=xi(i);
u(i)=1;

else
x(i)=xi(i);
u(i)=0;
end
end
plot(x,u,'b');
xlabel('x');           
  ylabel('u(x,t)');
hold on
t=0.25;
for i=1:2000
    if xi(i)< t
    x(i)=xi(i);
    u2(i)=x(i)/t;


    elseif 1+t/2<xi(i)
     x(i)=xi(i);
     u2(i)=0;
    elseif t/2 <xi(i)< 1+t/2
     x(i)=xi(i);
     u2(i)=1;

end
end
plot(x,u2,'r');
xlabel('x');           
  ylabel('u(x,t)');
hold on  xlabel('x');
t=0.5;
for i=1:2000
    if  xi(i)< t
    x(i)=xi(i);
    u3(i)=x(i)/t;


    elseif 1+t/2<xi(i)
     x(i)=xi(i);
     u3(i)=0;
    elseif t/2 <xi(i)< 1+t/2
     x(i)=xi(i);
     u3(i)=1;

end
end
plot(x,u3,'g');
xlabel('x');           
  ylabel('u(x,t)');
axis([0,pi,0,1.25])
hold on  xlabel('x');           
  ylabel('u(x,t)');
  legend('T=0','T=0.25','T=0.5')
  title('1D Burgers Analytic solution at time periods of 0.25');
  hold on