% Numerical experiments
% The Numerical experiments of this paper:
%% The Numerical experiment in A.1.1 experiment on $C_2^{25}$}:
Demo_Random_example
% Table is the result
% rows of Table is: [Lower bound ,time costs, degree of approx polynomial ,number of basis];
%% The Numerical experiment in A.1.2 experiment on $C_3^{15}$:
Demo_Random_example2
% Table is the result
% rows of Table is: [Lower bound ,time costs, degree of approx polynomial ,number of basis];
%% The Numerical experiment  A.2.1 experiment on unweighted MAX-2SAT problems and experiment on weighted MAX-3SAT problems.
% Set the cnf/ wcnf File into the "dir" directory, 
% set the d=degree of approx Polynomial ;
% number of basis= sparsity* |supp(f)|;
% here is an example of 
%d = 2, k = 1.5 |supp(f)|,?l = 0, m = sum of  weights in A.2.2
d=2;
sparsity=1.5;
Running_Dir;
% Table is the result
% rows of Table is: [Lower bound ,time costs, degree of approx polynomial ,number of basis];

%%  The Numerical experiment in A.3.2 experiment on benchmark problems.
T=[];
File='1.wcnf';
tol=1e-7;
rho_fun=@(rho,N,i)rho.^N;
rho_fun2=@(rho,N,i)2.^(-(i-1));
[lb,Q,y,x,f,Index]=HLM08(File,tol,{rho_fun,rho_fun2},'p');
y_rho_i_N=y(1);
y_2_i=y(2);
k=size(Index,1);
[Q,f,Index,y_Gram,y_Moment,lb]=OurRounding(File,k);
T=[y_rho_i_N,y_2_i,y_Gram,y_Moment];
% T is the result
% Rows of Table is:  the rounding results of [rho_i^N, 2^{-(i-1)}, Gram,  Moment]
