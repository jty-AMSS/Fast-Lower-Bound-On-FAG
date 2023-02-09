function [blk,At,C,b]=CVX2SDPNALPlus(Ap,Aq)

%Ap*x(:)=Aq
Np=size(Ap,2);
Np=sqrt(Np);
N=Np;
I=(1:N^2);I= reshape(I,[N,N]);
I=I';
Tril_m1=tril(I,-1);Tril_m1=Tril_m1';
Tril_m1=Tril_m1(:);
Tril_m1=Tril_m1(Tril_m1>0);
Tril_0=tril(I);Tril_0=Tril_0';
Tril_0=Tril_0(:);
Tril_0=Tril_0(Tril_0>0);

Ap(:,Tril_m1)=sqrt(2)*Ap(:,Tril_m1);
Ap=Ap(:,Tril_0);
Ap=Ap.';
ObjFun=zeros(Np);
ObjFun(1)=-1;
C={ObjFun};
At={Ap};
b=Aq;
blk={'s',[Np]};
%  [obj,X,s,y,S,Z,ybar,v,info,runhist]= sdpnalplus(blk,At,C,b,[],[],[],[],[],opti);


end