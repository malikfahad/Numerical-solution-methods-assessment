clear all;
clc;
format long
%Inputs: Interior points
m=64;
%Total points
Mt=m+2;
%h
h=1/(m+1);
%Grid
x=(0:h:1);
y=(0:h:1);

tol=1d-6;
%Intialization of the Matrices
N=Mt*Mt;              %N Unknowns
A=sparse(zeros(N,N)); %Matrix A written in block-diagonal form
F=sparse(zeros(N,1)); %Vector F is the result. The boundary conditions are absorbed into F.
%Interior Points
for i=2:Mt-1 %Running Loop over x, the first and last grid points (Boundary condtions) are omitted
    for j=2:Mt-1 %Running Loop over y, the first and last grid points (Boundary condtions) are omitted
        n=i+(j-1)*Mt; %Enumeration of the unknown terms in Vector V
        A(n,n)=-4; %Diagonal 
        A(n,n-1)=1; %sub-diagonal
        A(n,n+1)=1; %super-diagonal
        A(n,n-Mt)=1;%far off sub-diagonal
        A(n,n+Mt)=1;%far off super-diagonal
        F(n,1)=-2*x(j);
    end
end


% Boundary Points: %BC: u(0,y)=0
i=1;
for j=1:Mt
    n=i+(j-1)*Mt;
    A(n,n)=1;
    F(n,1)=0;
end
%BC: u(0,y)=0
i=Mt;
for j=1:Mt
    n=i+(j-1)*Mt;
    A(n,n)=1;
    F(n,1)=0;
end
%BC: u(x,0)=0
j=1;
for i=1:Mt
    n=i+(j-1)*Mt;
    A(n,n)=1;
    F(n,1)=0;
end
%BC: u(x,1)=0
j=Mt;
for i=1:Mt
    n=i+(j-1)*Mt;
    A(n,n)=1;
    F(n,1)=0;
end
maxit=m*m;

D=diag(diag(A));
U=triu(A)-diag(diag(A));
L=tril(A)-diag(diag(A));


method = menu('Choose a Gradient Method','Steepest Descent','Conjugate Gradient');
if method==1
    technique = menu('Choose a technique','No Preconditioning','Preconditioning');
    if technique == 2
        preconditioning = menu('Choose a preconditioning','Jacobi iteration', 'Symmetric Gauss–Siedel','Symmetric SOR');
        if preconditioning == 1
            M1=D;
            M2=speye(size(A));
            M=M1*M2;
            [k, u,e] = pcg3(A,F,M,method,tol,maxit);
        elseif preconditioning == 2
            M1=(D+L);
            M2=inv(D)*(D+U);
            M=M1*M2;
            [k, u,e] = pcg3(A,F,M,method,tol,maxit);
        else
            M=A;
            [k, u,e] = pcg3(A,F,M,method,tol,maxit);
        end
    else
        M=speye(size(A));
        [k, u,e] = pcg3(A,F,M,method,tol,maxit);
    end
 
else 
    technique = menu('Choose a technique','No Preconditioning','Preconditioning');
    if technique == 2
        preconditioning = menu('Choose a preconditioning','Jacobi iteration', 'Symmetric Gauss–Siedel','Symmetric SOR');
        if preconditioning == 1
            M1=D;
            M2=speye(size(A));
            M=M1*M2;
            [k, u,e] = pcg3(A,F,M,method,tol,maxit);
        elseif preconditioning == 2
            M1=(D+L);
            M2=inv(D)*(D+U);
            M=M1*M2;
            [k, u,e] = pcg3(A,F,M,method,tol,maxit);
        else
            M=A;
            [k, u,e] = pcg3(A,F,M,method,tol,maxit);
        end
    else
        M=speye(size(A));
        [k, u,e] = pcg3(A,F,M,method,tol,maxit);
    end
end


for i=1:Mt
    for j=1:Mt
        n=i+(j-1)*Mt;
        V(i,j)=u(n);
    end
end


figure
mesh(0:h:1,0:h:1,V)
title(['Numerical Solution for m = ' num2str(m)'']);
xlabel('x')
ylabel('y')

fprintf('Number of Iterations: %d \n',k)
fprintf('CPU Time: %d \n',e)

    
