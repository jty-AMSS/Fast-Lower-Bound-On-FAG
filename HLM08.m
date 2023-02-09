function [lb,Q,y,x,f,Index]=HLM08(File,tol,rho_fun,BasisStr)
% 复现了HLM08的rounding 程序：
% input：
% File：cnf/wcnf path
% tol: tol of eign=0 (e.g. 1e-7);
% rho_fun function handle, which is skewed factor of [HLM08]
% vector computation method， rho=(1-\mu_i+mu_1)，N=corank of solution of SDP，i=sort of eignvalue，e.g：rho_fun={@(rho,N,i)rho.^N,@(rho,N,i)2.^(-(i-1))}
% % out:
% lb: SDP lower bound
% Q：SDP solution
% y: min (f) by rounding
% x:argmin(f) by rounding=

if nargin==3
    BasisStr='p';
end
[Index,f]=Basis08(File,BasisStr);
m=size(Index,1);
%% SDP solve
[A,B,~,~,Suppf_lb]=GenerateConidition(f,Index);
cvx_begin
variable  Q(m,m) semidefinite
minimize(trace(Q) )
A*Q(:)==B;
cvx_end
n=m;
F=@(x)trace(0.5*reshape(x,[n,n])+0.5*reshape(x,[n,n])')-n*lambda_min(0.5*reshape(x,[n,n])+0.5*reshape(x,[n,n])')+norm(A*x(:)-B,1);
lb= f(zeros(1,length(f.n)))-F(Q)-Suppf_lb;

%% Rounding
Deg_1_Index=[];
Deg_2_Index=[];
for i=1:m
    t=Index(i,:);t=t+0;
    if sum(t)==0
        Deg_0_Index=i;
    end
    if sum(t)==1
        Deg_1_Index(find(t))=i;
    end
    if sum(t)>1
        Deg_2_Index=[Deg_2_Index;i];
    end
end
Q=full(Q);
[x,mu] = eig(Q,'vector');
V_Index=(abs(mu)<tol);
V=x(:,V_Index);
B=V(Deg_2_Index,:);
BB=B'*B;
[U,mu] = eig(BB,'vector');
mu1=min(mu);
rho=1-(mu-mu1);
N=length(BB);
[~,i]=sort(mu);
for index_rho_fun=1:length(rho_fun)
    rho_fun0=rho_fun{index_rho_fun};
    w=rho_fun0(rho,N,i);
    w=w(:);
    w0=randn(N,1);
    w0=w0./norm(w0);
    w=w.*w0;
    lambda=U*w;
    x=V*lambda(:);
    signx=sign(x(Deg_0_Index));
    x=x(Deg_1_Index);
    x=sign(x)*signx;
    x(x==1)=0;
    x(x==-1)=1;
    y(index_rho_fun)=f{x};
end
% V=CZifft(f);
% minv=min(V(:));
end




function [Index,f,w]=Basis08(File,BasisStr)
% BasisStr='GW' ,'p','ap' 't' 'pt'
if nargin<2
    BasisStr='p';
end
[f,w,D]=DICMS2function(File);
Index=PolySOSWithAlpha(f,D,BasisStr);
end
function subs= PolySOSWithAlpha(f,D,BasisStr)
n=f.n;n=length(n);
subs=[];
subs=[subs;char(zeros(1,n))];
subs=[subs;char(eye(n))];
%GW end
if (strcmp(BasisStr,'p')||strcmp(BasisStr,'pt'))
    MpIndex=[];
    Mp=[];
    for i=1:length(D)
        t=D{i};
        if (length(t)>=2)
            t=abs(t);
            t=sort(t);
            MpIndex=[MpIndex;nchoosek(t,2)];
        end
    end
    MpIndex=unique(MpIndex,'rows');
    Mp=zeros(size(MpIndex,1),n);
    for i=1:size(MpIndex,1)
        Mp(i,MpIndex(i,:))=1;
    end
    Mp=char(Mp);
    subs=[subs;Mp];
end
%M_p End
if (strcmp(BasisStr,'t')||strcmp(BasisStr,'pt'))
    MtIndex=[];
    Mt=[];
    for i=1:length(D)
        t=D{i};
        if (length(t)>=3)
            t=abs(t);
            t=sort(t);
            MtIndex=[MtIndex;nchoosek(t,3)];
        end
    end
    MtIndex=unique(MtIndex,'rows');
    Mt=zeros(size(MtIndex,1),n);
    for i=1:size(MtIndex,1)
        Mt(i,MtIndex(i,:))=1;
    end
    Mt=char(Mt);
    subs=[subs;Mt];
end
%M_t End
if (strcmp(BasisStr,'ap'))
    MapIndex=nchoosek(1:(f.n),2);
    for i=1:size(MapIndex,1)
        M_ap(i,MapIndex(i,:))=1;
    end
    for i=1:size(MapIndex,1)
        M_ap(i,MapIndex(i,:))=1;
    end
    M_ap=char(M_ap);
    subs=[subs;M_ap];
end
%M_ap End
if (strcmp(BasisStr,'findf'))
    subs=[subs;find(f)];
end
subs=unique(subs,'rows');
end

function z=lambda_min(Q)
z=eig(Q);
z=min(z);
end


