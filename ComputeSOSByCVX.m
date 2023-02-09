function [P,A,b]=ComputeSOSByCVX(f,Index)
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
cvx_begin
variable  P(m,m) hermitian semidefinite
minimize(norm(A*P(:)-b,1)+trace(P) )%-length(P)*lambda_min(P))
%A*P(:)==b;
cvx_end

end
