
%% on Z_2^25:
d=2;
n=25;
degree_bound=3;
num_terms=450;
Table=[];
F={};M=[];
for i=1:10
    [f,m]=Random_Group_function(ones(1,n)+1,num_terms,degree_bound);
    F{i}=f;
    M(i)=m;
    sparsity=3*length(f);
    tic;
    [Q,~,Index,A,B,~,Suppf_lb,X0]=FastLowerBound_Poly(f,2,(0:m)',sparsity);
    lb_l=norm(A*Q(:)-B,1)+trace(Q)- length(Q)*lambda_min(Q);
    f_const=f((zeros(1,length(f.n))));
    lb_l=f_const-lb_l;
    lb_l=lb_l-Suppf_lb;
    time=toc;
    tline=[lb_l,time,d,size(Index,1)];
    Table=[Table;tline];
%     disp(vpa(tline,5))
%     txtname=['n',num2str(n),'t',num2str(num_terms),'d',num2str(degree_bound),'_',num2str(i),'.jl']
%     fid=fopen(['random_example\' txtname],'w');
%     fprintf(fid,'using TSSOS\n');
%     fprintf(fid,'using DynamicPolynomials\n');
%     fprintf(fid,['@polyvar x[1:', num2str(n),'] \n']);
%     str=CZ2TSSOS(f);
%     fprintf(fid,['f=' str,'\n']);
%     fprintf(fid,['d=2;import Dates;x1=time();opt,sol,data = tssos_first([f],nb=', num2str(n) ',x, d, TS="MD");x2=time();opt0,sol,data = tssos_higher!(data, TS="MD");x3=time();T21=x2-x1;T31=x3-x1;']);
%     fprintf(fid,['x1=time();opt,sol,data = cs_tssos_first([f],nb=', num2str(n) ',x, d, TS="MD");x2=time();opt1,sol,data = cs_tssos_higher!(data, TS="MD");x3=time();T21cs=x2-x1;T31cs=x3-x1;']);
%     fprintf(fid,['print([T21,T31,opt0,T21cs,T31cs,opt1])']);
%     fclose(fid);
end

