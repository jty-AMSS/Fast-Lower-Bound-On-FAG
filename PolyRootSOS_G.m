function [Q,Index,A,B,g,F,Suppf_lb]=PolyRootSOS_G(f,d,range,k)
%range: range of function f£¬class£º 1-d array
p=range;
Poly=Polyfit_inf_norm(p,d);
g=Poly(2)*f+Poly(1);
fd=f;
for i=3:length(Poly)
    fd=fd.*f;
    g=Poly(i)*fd+g;
end
[a,b]=find(g);
[~, sortb_Index]=sort(abs(b),'descend');
k=min(k,size(a,1));
% k=ceil(sqrt(length(f)));
Index= a(sortb_Index(1:k),:);
[A,B,~,~,Suppf_lb]=GenerateConidition(f,Index);
[X0,lb_afpc] = AFPC_BB_moment(size(Index,1),A,B,1e3,1e-4,1/4,1,'full');
disp(lb_afpc)
Q=ComputeSOSByGrid_Linear(X0,A,B,Suppf_lb);
F=@(x)trace(0.5*reshape(x,[n,n])+0.5*reshape(x,[n,n])')-n*lambda_min(0.5*reshape(x,[n,n])+0.5*reshape(x,[n,n])')+norm(A*x(:)-B,1);
end


function z=lambda_min(Q)
z=eig(Q);
z=min(z);
end

function Poly=Polyfit_inf_norm(x,d)
A=[x.^[0:d],-ones(length(x),1);-x.^[0:d],-ones(length(x),1)];
Poly=linprog([zeros(1,d+1),1],A,[x.^0.5;-x.^0.5]);
Poly=Poly(1:d+1);
end


