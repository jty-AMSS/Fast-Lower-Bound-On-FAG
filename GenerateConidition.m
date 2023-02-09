function [A,b,Eq,Suppf_Minus_Condition,Suppf_lb]=GenerateConidition(f,Index,ifconsider_complex)
if nargin==2
    ifconsider_complex=1;
end
if all(f.n==2)
    
    Eq=containers.Map('KeyType',  'char', 'ValueType', 'any');
    m=size(Index,1);
    N=f.n;
    for i=1:m
        for j=1:m
            t=mod(Index(i,:)-Index(j,:),N);
            if isKey(Eq,char(t))
                Eq(char(t))=[Eq(char(t)),m*(i-1)+j];
            else
                Eq(char(t))=[m*(i-1)+j];
            end
        end
    end
    
    suppf=find(f);
    G=keys(Eq);
    %%%
    Eq.remove(char(zeros(1,length(f.n))));
    %%%
    G=keys(Eq);
    Ax=[];Ay=[];b=zeros(length(G),1);
    for i=1:length(G)
        t=G{i};
        Ax(end+1:end+(length(Eq(t))))=i;
        Ay(end+1:end+(length(Eq(t))))=Eq(t);
        b(i)=f(t);
    end
    A=sparse(Ax,Ay,1,length(G),m*m);
    G{end+1}=char(zeros(1,length(f.n)));
    Suppf_Minus_Condition=setdiff(suppf,cell2mat(G(:))+0,'rows');
    Suppf_lb=0;
    for i=1:size(Suppf_Minus_Condition,1)
        Suppf_lb=Suppf_lb+abs(f(char(Suppf_Minus_Condition(i,:))));
    end
    return
end
if ~isreal(b)&&ifconsider_complex
    [A,b]=complex2real(A,b);
end

end