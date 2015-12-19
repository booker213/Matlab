clc 
clear
format long
%% Define Constants
q=0.99; alpha=0.4709; beta = 1;
H=1; Bbase= 1.178164343; Bt=0.585373798;


m=11;  n=m+1;
dz=1/m; dt=0.0001;
mu=dt/dz^2;

test=dz^2/(alpha*dz*Bbase^2+2*beta*Bbase^3)

%% Define Variables

z_steady=linspace(0,H,n);
b_left=zeros(n,1);
b_right=zeros(n,1);
b_highres=zeros(n,1);
z= linspace(0,H,m+1);
b=zeros(m+1,1);

%% Steady State


dz1=H/n;

b_steady=Bbase;

%RK4 scheme for steady state
for i=1:n
k1= dz1*Numhw1deriv(b_steady);
k2= dz1*Numhw1deriv(b_steady+0.5*k1);
k3= dz1*Numhw1deriv(b_steady+0.5*k2);
k4= dz1*Numhw1deriv(b_steady+k3);
b_steady= b_steady+(1/6)*(k1+2*k2+2*k3+k4);
b_left(i)=-b_steady/2;
b_right(i)=-b_left(i);
b_highres(i)=2*b_right(i);
end

error_rk= Bt-b_steady


figure
plot(b_left,z_steady,'b')
hold on
plot(b_right,z_steady,'b')
axis([-1 1 0 1])
title('Steady State Dyke Width')
xlabel('Dyke Width - b(z)')
ylabel('z')

%% Finite Difference Solution

%Set Up Initial Condition
b(1)=Bbase;
for i=2:m+1
b(i)=b_steady;
end
bnew=b;
for time=0:dt:2
% Finite Difference scheme for adjoint form
for i=2:m
bnew(i)=b(i)-mu*dz*alpha*(b(i)^3-b(i-1)^3) + beta*mu*((0.5*(b(i+1)+b(i)))^3 ...
*(b(i+1)-b(i)) + (0.5*(b(i)+b(i-1)))^3*(b(i-1)-b(i)));

end
b=bnew;
% If statements to output dimensional vectors at required times.
if time==0.05
b_005=b;
b_left005=-b/2;
b_right005=b/2;
end

if time==0.1
b_01=b;
b_left01=-b/2;
b_right01=b/2;
end

if time==0.2
b_02=b;
b_left02=-b/2;
b_right02=b/2;
end

if time==0.5
b_05=b;
b_left05=-b/2;
b_right05=b/2;
end

if time==1
b_1=b;
b_left1=-b/2;
b_right1=b/2;
end
if time==2
b_2=b;
b_left2=-b/2;
b_right2=b/2;
end

end

%Plot dimensional graphs
figure
%Plot left side first for legends
plot(b_left,z_steady,'b')
hold on
plot(b_left005,z,'r')
hold on
plot(b_left01,z,'g')
hold on
plot(b_left02,z,'k')
hold on
plot(b_left05,z,'c')
hold on
plot(b_left1,z,'m')
hold on
plot(b_left2,z,'y')
hold on
plot(b_right,z_steady,'b')
hold on
plot(b_right005,z,'r')
hold on
plot(b_right01,z,'g')
hold on
plot(b_right02,z,'k')
hold on
plot(b_right05,z,'c')
hold on
plot(b_right1,z,'m')
hold on
plot(b_right2,z,'y')
axis([-1 1 0 1])
title([num2str(m), ' point Dyke Width'] );
legend('Steady','time = 0.05','time = 0.1','time=0.2','time=0.5', ...
    'time=1','time=2')
xlabel('Dyke Width - b(z,t)')
ylabel('z')


%% Errors

error_L2 =((b_highres(1)-b(1))^2) ;
for i=2:n
error_L2 = error_L2 + ((b_highres(i)-b(i))^2);
end
error_L2 = sqrt((H/(2*m))*error_L2)

error_max = abs(b_highres(1)-b(1));
for i=2:n
   if abs(b_highres(i)-b(i))> error_max
      error_max = abs(b_highres(i)-b(i));
   end
end
error_max
   

%% Exact Solution
c=alpha; z_0=0;
