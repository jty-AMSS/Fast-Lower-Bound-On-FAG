function [lb_l,Q,f,X0,Index,x,y]=FSOSBulider(path,d,sparsity,opti)
if nargin==3
    opti=0;
end


[f,w,~]=DICMS2function(path);
if sparsity<3
    sparsity=sparsity-1e-10;
end
if abs(sparsity-round(sparsity))>0
    sparsity=length(f)*sparsity;
    sparsity=round(sparsity);
end
if sparsity==-1
    sparsity=length(f);
end


range=0:sum(w);range=range(:);
[Q,~,Index,A,B,~,Suppf_lb,X0]=FastLowerBound_Poly(f,d,range,sparsity);
lb_l=norm(A*Q(:)-B,1)+trace(Q)- length(Q)*lambda_min(Q);
f_const=f((zeros(1,length(f.n))));
lb_l=f_const-lb_l;
lb_l=lb_l-Suppf_lb;
end