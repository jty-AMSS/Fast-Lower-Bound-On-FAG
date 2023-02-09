function [A,b,Eqi,Suppf_Minus_Condition,Suppf_lb]=GenerateConidition_complex(f,Index)
Eqi=containers.Map('KeyType',  'char', 'ValueType', 'any');
Eqj=containers.Map('KeyType',  'char', 'ValueType', 'any');
m=size(Index,1);
N=f.n;
for i=1:m
    for j=1:m
        t=mod(Index(i,:)-Index(j,:),N);
        if isKey(Eqi,char(t))
            Eqj(char(t))=[Eqj(char(t)),i];
            Eqi(char(t))=[Eqi(char(t)),j];
        else
            Eqj(char(t))=[i];  Eqi(char(t))=[j]; %(j,i)->m*(i-1)+j
        end
    end
end

suppf=find(f);
G=keys(Eqi);
%%%
Eqi.remove(char(zeros(1,length(f.n))));

Eqj.remove(char(zeros(1,length(f.n))));
%%%
G=keys(Eqi);
Ax=[];Ay=[];b=zeros(2*length(G),1);
Axm=[];Aym=[];
G_m=length(G);
for i=1:length(G)
    t=G{i};
    Ax(end+1:end+(length(Eqi(t))))=i;
    Ay(end+1:end+(length(Eqi(t))))=Eqi(t)+2*m*(Eqj(t)-1); %(i,j)_2m
    Ax(end+1:end+(length(Eqi(t))))=i;
    Ay(end+1:end+(length(Eqi(t))))=m+Eqi(t)+2*m*(Eqj(t)-1+m);%(i+m,j+m)
    
    Ax(end+1:end+(length(Eqi(t))))=i;
    Ay(end+1:end+(length(Eqi(t))))=m+Eqj(t)+2*m*(Eqi(t)-1+m);%(j+m,i+m)
    Ax(end+1:end+(length(Eqi(t))))=i;
    Ay(end+1:end+(length(Eqi(t))))=Eqj(t)+2*m*(Eqi(t)-1);%(j,i)
    
% % % % % % % % % % %     
    Axm(end+1:end+(length(Eqi(t))))=length(G)+i;
    Aym(end+1:end+(length(Eqi(t))))=m+Eqi(t)+2*m*(Eqj(t)-1);%(i+m,j)
    Ax(end+1:end+(length(Eqi(t))))=length(G)+i;
    Ay(end+1:end+(length(Eqi(t))))=Eqi(t)+2*m*(Eqj(t)-1+m);%(i,m+j)
    
    Axm(end+1:end+(length(Eqi(t))))=length(G)+i;
    Aym(end+1:end+(length(Eqi(t))))=Eqj(t)+2*m*(m+Eqi(t)-1);%(j,i+m)
    Ax(end+1:end+(length(Eqi(t))))=length(G)+i;
    Ay(end+1:end+(length(Eqi(t))))=m+Eqj(t)+2*m*(Eqi(t)-1);%(m+j,i)
    b(i)=real(f(t));
    b(G_m+i)=imag(f(t));
end
Ap=sparse(Ax,Ay,1/2,2*length(G),4*m*m);
Am=sparse(Axm,Aym,-1/2,2*length(G),4*m*m);
A=Ap+Am;
G{end+1}=char(zeros(1,length(f.n)));
Suppf_Minus_Condition=setdiff(suppf,cell2mat(G(:))+0,'rows');
Suppf_lb=0;
for i=1:size(Suppf_Minus_Condition,1)
    Suppf_lb=Suppf_lb+abs(f(char(Suppf_Minus_Condition(i,:))));
end
end