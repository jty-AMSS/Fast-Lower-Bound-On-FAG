function [c,x]=RoundingByMoment(f,Index)
Index=union(Index,eye(size(Index,2)),'rows');
if isreal(f)
    [A,B,~,~,Suppf_lb]=GenerateConidition(f,Index);
else
    [A,B,~,~,Suppf_lb]= GenerateConidition_complex(f,Index);
end
[obj,X,s,y,S,Z,ybar,v,info,runhist]=ComputeSOSBySDPNAL2(A,B);
S=S{1};
m=size(Index,1);Deg_2_Index=[];
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
[u,v]=eig(S,'vector');
v(1:end-1)=0;
S0=u*diag(v)*u';
S0=S0+S0';
S0=S0/2;
x=S0(1,Deg_1_Index);
x=sign(x)*sign(S0(1,Deg_0_Index));
x(x==1)=0;x(x==-1)=1;
c=f{x};
end


function  [obj,X,s,y,S,Z,ybar,v,info,runhist]=ComputeSOSBySDPNAL2(A,b)
m=size(A,2);
m=sqrt(m);
opti=SDPNALplus_parameters;
opti.maxtime=1e6;
ADMmaxiter=max(round(m),300);
opti.ADMmaxiter=ADMmaxiter;
opti.maxiter=1;
[blk,At,C,b]=CVX2SDPNALPlus(A,b);
C{1}=eye(m);
if isreal(b)
    [obj,X,s,y,S,Z,ybar,v,info,runhist]= sdpnalplus(blk,At,C,b,[],[],[],[],[],opti);
else
    [bblk,AAt,CC,bb,iscmp] = convertcmpsdp(blk,At,C,b);
    [obj,X,s,y,S,Z,ybar,v,info,runhist]= sdpnalplus(bblk,AAt,CC,bb,[],[],[],[],[],opti);
end
end