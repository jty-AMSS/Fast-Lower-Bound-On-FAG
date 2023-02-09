function [yN,yi,ys,ym,f]=random_MAX_SAT(n,d,FileIndex)
if nargin==2
    FileIndex=[];
end
% Generate Rounding_example\randomwcnf .wcnf 
% n=number of vars
% number of clauses= d* n
% yN: rho_i^N   method in [HLM08]
% yN: 2^(-(i-1)) method in [HLM08]
% ys: corank1 Gram-matrix based
% ym £ºraml one approx moment matrix based

m=d*n;
File=['Rounding_example\randomwcnf' num2str(FileIndex) '.wcnf'];
% File='randomwcnf.wcnf';
fid = fopen(File,'w+');
fprintf(fid,'c Weighted CNF');
fprintf(fid,'\n');
fprintf(fid,'c from my generator');
fprintf(fid,'\n');
fprintf(fid,'p wcnf  ');
fprintf(fid,num2str(n));
fprintf(fid,'  ');
fprintf(fid,num2str(m));
for i=1:m
    fprintf(fid,'\n');
    t= randperm(n,2);
    y=randn(1,2);y=sign(y);
    y=y.*t;
    y=[1,y,0];
    fprintf(fid,num2str(y));
end
fclose(fid);
% end



tol=1e-7;
rho_fun=@(rho,N,i)rho.^N;
rho_fun2=@(rho,N,i)2.^(-(i-1));
[lb,Q,y,x,f,Index]=HLM08(File,tol,{rho_fun,rho_fun2},'p');
yN=y(1);
yi=y(2);
% [lb,Q,yi,x,f,Index]=HLM08(File,tol,rho_fun2,'p');
k=size(Index,1);
[Q,f,Index,ys,ym,lb]=OurRounding(File,k);
end