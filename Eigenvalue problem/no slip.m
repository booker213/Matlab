%Basic code to determine the largest eigenvalue of f''+ sf = 0, f(0)=f(1)=0
%Discretise with finite differences on a grid of N points x=[x1,...,xN]
%Generalised eigenvalue problem of the form Ax=sBx, where A,B are NxN matrices

function [C,I] = free(Ra, N, a, Pr)

%%Initialise block matrices
A=zeros(3*N,3*N);
B=zeros(3*N,3*N);

%%boundary conditions (here, two BCs so they fill the first two rows).
A(1,1)=1;
A(2,N)=1;
A(3,1)=-(N-1)*3/2;
A(3,2)=(N+1)*4/2;
A(3,3)=-(N-1)/2;
A(4,N-2)=(N-1)/2;
A(4,N-1)=-(N-1)*4/2;
A(4,N)=(N-1)*3/2;
A(5,2*N+1)=1;
A(6,3*N)=1;

%%inner matrix elements of A and B

for j=1:N-2
A(6+j,j)= -a^2*(N-1)^2;    
A(6+j,j+1)= a^2*(2*(N-1)^2 + a^2);
A(6+j,j+2)=-a^2*(N-1)^2; 
A(6+j,N+j)= (N-1)^2;
A(6+j,N+j+1)= -(2*(N-1)^2 + a^2);
A(6+j,N+j+2)= (N-1)^2;
A(6+j,2*N+j+1)= -a^2*Ra;
A(6+(N-2)+j,j+1)= 1;
A(6+(N-2)+j,2*N+j)= (N-1)^2;
A(6+(N-2)+j,2*N+j+1)= -(2*(N-1)^2 + a^2);
A(6+(N-2)+j,2*N+j+2)= (N-1)^2;
A(6+2*(N-2)+j,j)= (N-1)^2;
A(6+2*(N-2)+j,j+1)= -2*(N-1)^2;
A(6+2*(N-2)+j,j+2)= (N-1)^2;
A(6+2*(N-2)+j,N+j+1)= -1;

B(6+j,j)= (Pr^(-1))*(N-1)^2;
B(6+j,j+1)= -(Pr^(-1))*(2*(N-1)^2 + a^2);
B(6+j,j+2)= (Pr^(-1))*(N-1)^2;
B(6+(N-2)+j,2*N+j+1)= 1;
end

%if several block matrices to form A and B,
%concatenate them with a command like 
%A=vertcat([A11 A21 A31],[A12 A22 A32],[A13 A23 A33]);

%calculate eigenvalues
%e=eig(A,B,'qz')
%here, we use the function qz, which allows an intermediate step, where we eliminate spurious eigenvalues

format long e
[AA,BB,Q,Z,V,W]=qz(A,B);
alpha=diag(AA);
beta=diag(BB);

j=1;
for i=1:length(alpha)
if beta(i)>0 | beta(i)<0
    e(j)=alpha(i)/beta(i);
    F(j,:)=V(:,i);
    j=j+1;
end
end

%display(e)

%max eigenvalue
[C,I]=max(e)

%associated eigenvector
ev=F(I,:);hoo
plot((1:N)/(N-1), ev(1:N))
xlabel 'x'
ylabel 'eigenfunction'
xlim([0 1])

clear A AA B BB Q Z V W

