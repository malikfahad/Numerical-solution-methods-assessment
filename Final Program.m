clc
clear;
clear all;
%% Initialisation of parameters
% Number of grid points to record number of iteration needed
n1 = 10; 
n2 = 20;
n3 = 40;

init_point = 0; % Initial point
end_point = 1; % End point
range = end_point-init_point; % Interval length
tol = 1e-5; % Tolerance level
maxit = 4000; % Maximum number of iterations
i = 1; % Index 

%% For-loop to test on gradient(mg) and conjugate gradient(cg) method
for n = 10:n3
    if n == n1 || n == n2 || n == n3
        % Calculation of h and x,y, depending on endpt
        h = range/(n+1);
        x = init_point:h:end_point;
        y = init_point:h:end_point;
        
        [xnew,ynew] = meshgrid(x,y);
        
        f = -2*xnew;
        
        gb = 0*x + 0*y; % u(x,0)
        gt = 0*x + 0*y; % u(x,endpoint)
        gl = 0*x + 0*y; % u(0,y)
        gr = 0*x + 0*y; % u(endpoint,y)
        
        A = Form_A(n);
        
        % Form_B.m forms the right hand side of the problem, B
        b = Form_B(n,h,f,gb,gt,gl,gr);
        
        % Call the DF sine FTransform routine for speed:
        b_T = FTransform(n,b)';
        [B_T_mod] = mat_mod(n,b_T);
        
        I = speye(n^2);
        
        Diag = diag(diag(A));
        
        C1 = tril(A);
        C2 = C1';
        
        weight = 2 - (2*pi/n);
        C3=((1-weight)/weight)*Diag+C1;
        C4=C3';
        
        tic;
        % CG with no precondition
        [u, iter(i,1),lg_iter(i,1)]=pcg_solver(A,B_T_mod,I,I,tol,maxit);
        time(i,1)=toc;
        lg_time(i,1)=log(time(i,1));
        
        tic;
        % PCG with Jacobi
        [u, iter(i,2),lg_iter(i,2)]=pcg_solver(A,B_T_mod,Diag,I,tol,maxit);
        time(i,2)=toc;
        lg_time(i,2)=log(time(i,2));
        
        tic;
        % PCG with symmetric Gauss-Seidel
        [u, iter(i,3),lg_iter(i,3)]=pcg_solver(A,B_T_mod,C1,C2,tol,maxit);
        time(i,3)=toc;
        lg_time(i,3)=log(time(i,3));
        
        tic;
        % PCG with symmetric SOR
        [u, iter(i,4),lg_iter(i,4)]=pcg_solver(A,B_T_mod,C3,C4,tol,maxit);
        time(i,4)=toc;
        lg_time(i,4)=log(time(i,4));
        
        tic;
        % MG with no precondition
        [u, iter(i,5),lg_iter(i,5)]=pmg_solver(A,B_T_mod,I,I,tol,maxit);
        time(i,5)=toc;
        lg_time(i,5)=log(time(i,5));
        
        tic;
        % PMG with Jacobi
        [u, iter(i,6),lg_iter(i,6)]=pmg_solver(A,B_T_mod,Diag,I,tol,maxit);
        time(i,6)=toc;
        lg_time(i,6)=log(time(i,6));
        
        tic;
        % PMG with symmetric Gauss-Seidel
        [u, iter(i,7),lg_iter(i,7)]=pmg_solver(A,B_T_mod,C1,C2,tol,maxit);
        time(i,7)=toc;
        lg_time(i,7)=log(time(i,7));
        
        tic;
        % PMG with symmetric SOR
        [u, iter(i,8),lg_iter(i,8)]=pmg_solver(A,B_T_mod,C3,C4,tol,maxit);
        time(i,8)=toc;
        lg_time(i,8)=log(time(i,8));
        
        i = i + 1;
    end
   
end

%% Graphical Plotting
x_map = [n1,n2,n3];

% Plot of Number of Iterations against N
figure(1);
plot(x_map,iter,'LineWidth',2)
title('Plot of Number of Iterations against N')
xlabel('N');
ylabel('Number of Iterations');
legend('CG','JCG','GSCG','SORCG','MG','JMG','GSMG','SORMG');

% Plot of ln(Number of Iterations) against N
figure(2);
plot(x_map,lg_iter,'LineWidth',2)
title('Plot of ln(Number of Iterations) against N')
xlabel('N');
ylabel('ln(Number of Iterations)');
legend('CG','JCG','GSCG','SORCG','MG','JMG','GSMG','SORMG');

% Plot of Elapsed Time against N
figure(3);
plot(x_map,time,'LineWidth',2)
title('Plot of Elapsed Time against N')
xlabel('N');
ylabel('Elapsed Time');
legend('CG','JCG','GSCG','SORCG','MG','JMG','GSMG','SORMG');

% Plot of ln(Elapsed Time) against N
figure(4);
plot(x_map,lg_time,'LineWidth',2)
title('Plot of ln(Elapsed Time) against N')
xlabel('N');
ylabel('lg(Elapsed Time)');
legend('CG','JCG','GSCG','SORCG','MG','JMG','GSMG','SORMG');











