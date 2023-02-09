function [Q,f,Index,ys,ym,lb,xs,xm]=OurRounding(File,num_basis)

[f,w]=DICMS2function(File);
d=2;
n=length(f.n);

range=0:sum(w);
p=range(:);
Poly=Polyfit_inf_norm(p,d);
g=Poly(2)*f+Poly(1);
fd=f;
for k=3:length(Poly)
    fd=fd.*f;
    g=Poly(k)*fd+g;
end
ADMmaxiter=max(round(num_basis));
for i=1:n
    t=zeros(n,1);
    t(i)=1;
    t=t(:);t=t';
    g(t)=inf;
end
t=zeros(n,1); t=t(:);t=t';
g(t)=inf;
[a,b]=find(g);
[~, sortb_Index]=sort(abs(b),'descend');

Index= a(sortb_Index(1:num_basis),:);
Index=union(Index,eye(size(Index,2)),'rows');
[A,b,~,~,Suppf_lb]=GenerateConidition(f,Index);
[blk,At,C,b]=CVX2SDPNALPlus(A,b);
Np=size(Index,1);
n=Np;
ObjFun=eye(Np);
C={ObjFun};
blk={'s',[Np]};
opti=SDPNALplus_parameters;
% opti.maxtime=maxtime;
opti.maxiter=1;
opti.ADMmaxiter=ADMmaxiter;
[obj,X,s,y,S,Z,ybar,v,info,runhist]= sdpnalplus(blk,At,C,b,[],[],[],[],[],opti);
P1=X{1};
disp(File)
disp('d')
n=length(P1);
x=P1(:);
lb(i)=trace(0.5*reshape(x,[n,n])+0.5*reshape(x,[n,n])')-n*lambda_min(0.5*reshape(x,[n,n])+0.5*reshape(x,[n,n])')+norm(A*x(:)-b,1);
[ys,xs]=Rounding(P1,Index,f);
[ym,xm]=RoundingByMoment(f,Index);
Q=P1;
end





function z=lambda_min(Q)
z=eig(Q);
z=min(z);
end

function Poly=Polyfit_inf_norm(x,d)
A=[x.^[0:d],-ones(length(x),1);-x.^[0:d],-ones(length(x),1)];
Poly=linprog([zeros(1,d+1),1],A,[x.^0.5;-x.^0.5]);
Poly=Poly(1:d+1);
if isnan(Poly(1))&&d==1
    Poly(1)=0;
    Poly(2)=1;
end
end
