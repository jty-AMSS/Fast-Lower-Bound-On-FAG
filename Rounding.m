function [y,x,Ind]=Rounding(Q,Index,f)
[a,mineig]= eig(Q,'vector');[~,k]=min(mineig);
v=a(:,k);
x=zeros(size(Index,2),1);
Ind=x;
Findflag(length(x),1)=0;
signFlag=1;
for i=1:size(Index,1)
    t=Index(i,:);
    if(sum(t)==0)
        signFlag=v(i);
    end
    if(sum(t)==1)
        t=find(t);
        Findflag(t)=1;
        x(t)=v(i);
        Ind(t)=i;
    end
end
NoValue=[];
if (sum(Findflag)<length(x))
    warning('not find:')
    disp(find(Findflag==0))
    for i=find(Findflag~=1)
        x(i)=-1;
    end
   NoValue=find(Findflag==0);
end
nNoValue=length(NoValue);
x0=sign(signFlag)*x;
x0=sign(x0);
x0(x0==1)=0;x0(x0==-1)=1;

y0=PointValue(f,x0);
for i=0:(2^nNoValue-1)
    t=dec2bin(i,nNoValue);
    t=t-'0';
    t(t==0)=-1;
    x0=sign(signFlag)*x;
    x0=sign(x0);
    x0(NoValue)=t;
    y=f{x0};
    if y<y0
        x=x0;
        y0=y;
    end
end
x=x0;
y=y0;
end