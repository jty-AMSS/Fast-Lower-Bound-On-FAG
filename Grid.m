function [g,A,b]=Grid(Q,A,b,n)
%min tr(Q)+||AQ(:)-b)||_1-len(Q)*min_eig(Q)
% [a,~]=eig(Q);
if nargin==3
n=length(Q);
end
y=A*Q(:)-b;
y2=A'*soft_sign(y,Q);
[a,mineig]= eig(Q,'vector');[~,k]=min(mineig);
tol=1e-16;
u=a(:,abs(mineig-mineig(k))<tol);
mineig=mineig(abs(mineig-mineig(k))<tol);
mineig=1-abs(mineig-mineig(k));
u=u*abs(mineig(:)+eps);
u=u/norm(u);
G1=-n*(u*u');

n=length(Q);
G2=reshape(y2,[n,n]);
G3=eye(n);
% G3(1)=0;
g=G1+G2+G3;
g=-g;
% g=-g;

end

function z=soft_sign(y,Q)
tol=1/length(Q)^2;
y1=y;
y(abs(y)>tol)=0;
y1=y1-y;
y1=y1/tol;
z=sign(y1)+y;
end

function z=hard_sign(y,Q)
tol=1/length(Q)^2;
y(abs(y)<tol)=0;
z=sign(y);
end