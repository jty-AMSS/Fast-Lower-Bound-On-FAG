function [AAp,bb]=complex2real(Ap,b)

%
% SDPNAL+ does not do complex SDP natively. To convert to real SDPs
% there are two approaches;
%   X >= 0 <==> [ real(X), -imag(X) ; imag(X), real(X) ] >= 0
%   X >= 0 <==> exists [ Y1, Y2^T ; Y2, Y3 ] >= 0 s.t.
%                        Y1 + Y3 == real(X), Y2 - Y2^T == imag(X)
% For primal standard form, this second form is the best choice,
% and what we end up using here.
%


br=real(b);
bi=imag(b);
bb=[br;bi];
m=size(Ap,2);
m=sqrt(m);
AAp=sparse(2*size(Ap,1),m^2*4);
for k=1:size(Ap,1)
    t=find(Ap(k,:));
    t=t(:);
    i=mod(t,m);
    i(i==0)=m;
    j=(t-i)/m+1;
    AAp(k,i+2*m*(j-1))=1;
    AAp(k,i+m+2*m*(j+m-1))=1;
    AAp(size(Ap,1)+k,i+m+2*m*(j-1))=1/2;
    AAp(size(Ap,1)+k,j+m+2*m*(i-1))=-1/2;
    AAp(size(Ap,1)+k,j+2*m*(i+m-1))=1/2;
    AAp(size(Ap,1)+k,i+2*m*(j+m-1))=-1/2;
end

end