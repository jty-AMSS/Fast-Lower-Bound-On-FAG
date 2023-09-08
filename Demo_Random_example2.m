
%% on Z_3^15:
d=2;%%%%%%%%%%%%%%%%%%%%%%
n=15;%%%%%%%%%%%%%%%%%%%%%%%%
degree_bound=2;%%%%%%%%%%%%%%%%%%%%%%%%5
num_terms=200;%%%%%%%%%%%%%%%%%%%%
Table=[];
F={};M=[];
for i=1:10
    [f,m]=Random_Group_function(ones(1,n)+2,num_terms,degree_bound);
    F{i}=f;
    M(i)=m;
    sparsity=round(3*length(f));%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    [Q,~,Index,A,B,~,Suppf_lb,X0]=FastLowerBound_Poly(f,d,(0:m)',sparsity);
    lb_l=norm(A*Q(:)-B,1)+trace(Q)- length(Q)*lambda_min(Q);
    f_const=f((zeros(1,length(f.n))));
    lb_l=f_const-lb_l;
    lb_l=lb_l-Suppf_lb;
    time=toc;
    tline=[lb_l,time,d,size(Index,1)];
    Table=[Table;tline];
%      %% CS_TSSOS
%     txtname=['n',num2str(n),'t',num2str(num_terms),'d',num2str(degree_bound),'_',num2str(i),'.jl']
%     fid=fopen(['random_example\' txtname],'w');
%     fprintf(fid,'using TSSOS\n');
%     fprintf(fid,'using DynamicPolynomials\n');
%     fprintf(fid,['@polyvar x[1:', num2str(2*n),'] \n']);
%     str=CZ2TSSOS(f);
%     fprintf(fid,['f=' str,'\n']);
%     for j=1:n
%         fprintf(fid,['g' num2str(j) '=2-x[',num2str(j),  ']^3- x[' num2str(j+n)  ']^3\n']);
%         fprintf(fid,['h' num2str(j) '=im*x[',num2str(j),  ']^3- im*x[' num2str(j+n)  ']^3\n']);
%     end
%     fprintf(fid,['pop=[f']);
%     for j=1:n
%         fprintf(fid,[',g' num2str(j), ',h'  num2str(j)]);
%     end
%     fprintf(fid,'];');
%     fprintf(fid,'d=4;import Dates;x1=time();\n ');
%     fprintf(fid,['n=' num2str(n) ';']);
%     fprintf(fid,['x1=time();opt,sol,data = cs_tssos_first(pop,x,n,d,numeq=', num2str(n*2) ',nb=n,TS="MD");x2=time();T21cs=x2-x1;print([opt,T21cs]);opt1,sol,data = cs_tssos_higher!(data, TS="MD");x3=time();T21cs=x2-x1;T31cs=x3-x1;']);
%     fprintf(fid,['print([opt,T21cs,opt1,T31cs])']);
%     fclose(fid);
end

