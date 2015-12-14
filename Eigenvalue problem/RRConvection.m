function [C,e,ev,I] = RRConvection(R,a,P,N)
%eigenvalue calculation: Given wavenumber, Rayleigh number, Prandtl number and number of
%grid points, calculates most unstable eigenvalue and associated
%eigenvector. 
    
%initialise block matrices Aij
A=zeros(3*N,3*N);
B=zeros(3*N,3*N);

%Boundary conditions
A(1,1)=-3*(N-1)/2;
A(1,2)=2*(N-1);
A(1,3)=-(N-1)/2;
A(2,N-2)=(N-1)/2;
A(2,N-1)=-2*(N-1); 
A(2,N)=3*(N-1)/2; 
A(3,1)= 1;
A(4,N)=1;
A(5,2*N+1)=1;   
A(6,3*N)=1;

%Inner matrix elements
for j=1:1:(N-2)
    A(6+j,j)=-a^2*(N-1)^2; 
    A(6+j,j+1)=a^2*(2*(N-1)^2+a^2); 
    A(6+j,j+2)=-a^2*(N-1)^2; 
    A(6+j,N+j)=(N-1)^2; 
    A(6+j,N+j+1)=-(2*(N-1)^2+a^2); 
    A(6+j,N+j+2)=(N-1)^2; 
    A(6+j,2*N+j+1)= -a^2*R;
    
    A(6+(N-2)+j,j+1)=1;
    A(6+(N-2)+j,2*N+j)=(N-1)^2; 
    A(6+(N-2)+j,2*N+j+1)=-(2*(N-1)^2+a^2); 
    A(6+(N-2)+j,2*N+j+2)=(N-1)^2; 
   
    A(6+2*(N-2)+j,j)=(N-1)^2; 
    A(6+2*(N-2)+j,j+1)=-2*(N-1)^2; 
    A(6+2*(N-2)+j,j+2)=(N-1)^2; 
    A(6+2*(N-2)+j,N+j+1)=-1;
    
    B(6+j,j)=(N-1)^2/P; 
    B(6+j,j+1)=-(2*(N-1)^2+a^2)/P; 
    B(6+j,j+2)=(N-1)^2/P;
    B(6+(N-2)+j,2*N+j+1)=1;
end

%QZ algorithm
[AA,BB,Q,Z,V,W]=qz(A,B,'real');
alpha=diag(AA);
beta=diag(BB);

%Eliminate spurious eigenvalues (N physical eigenvalues, but 3N eigenvalues
%in total
j=1;
for i=1:length(alpha)
if beta(i)>0 | beta(i)<0
    e(j)=alpha(i)/beta(i);
    F(j,:)=V(:,i);
    j=j+1;
end
end

%Select max eigenvalue and eigenvector (most unstable)
[C,I]=max(e)
ev=F(I,:);

%display(C);

%Plot eigenvector
plot(1:N, ev(1:N))
xlabel 'Grid point'
ylabel 'Eigenfunction'
%plot(1:N, F(:,1:N))

%Other eigenvalue function
%[v,e]=eig(A,B,'qz');
%plot(1:length(v),v);
%max(e(1:2*(N-2)))

clear A AA B BB Q Z V W
end

