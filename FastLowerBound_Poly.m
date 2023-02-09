function [Q,lb,Index,A,B,g,Suppf_lb,X0]=FastLowerBound_Poly(f,d,range,k)
%range: range of function £¬class£º 1-d array
% Isrounding=1 then add all linear terms to Index
p=range(:);
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
Index= a(sortb_Index(1:k),:);
if isreal(f)
    [A,B,~,~,Suppf_lb]=GenerateConidition(f,Index);
else
    [A,B,~,~,Suppf_lb]= GenerateConidition_complex(f,Index);
end
[X0,lb_nal]=ComputeSOSBySDPNAL(A,B);

g=Grid(X0,A,B);g=g+g';g=g/2;
% disp([lb_nal,norm(g(:))])
[Q]=ComputeSOSByGrid(X0,A,B,k);

n=length(Q);
F=@(x)trace(0.5*reshape(x,[n,n])+0.5*reshape(x,[n,n])')-n*lambda_min(0.5*reshape(x,[n,n])+0.5*reshape(x,[n,n])')+norm(A*x(:)-B,1);
lb= f(zeros(1,length(f.n)))-F(Q)-Suppf_lb;
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


