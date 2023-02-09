function [err]=CheckFSOS_complex(f,Index,Q)
n=length(Q);
m=n/2;
% Q=Q+eye(n)
% H=chol(Q);
Q=Q-min(eig(Q))*eye(n);
H=chol(sym(Q),'nocheck');
g=cell(1,n);
fsos=CZ(f.n);
for i=1:n
    try
    tempg=CZ(f.n);
    tempg2=CZ(f.n);
    for k=1:size(Index,1)
        tempg(char(Index(k,:)))=H(i,k)+1i*H(i,k+m);
    end
    g{i}=tempg;
    tempg2=tempg*conj(tempg);
    fsos=fsos+tempg2;
    catch
        1+1
    end
end
err=f-fsos;
end