function [Q,lb]=ComputeSOSBySDPNAL(A,b)
m=size(A,2);
m=sqrt(m);
opti=SDPNALplus_parameters;
opti.maxtime=1e6;
% opti.printlevel=2;
% opti.maxtime=2*m;
ADMmaxiter=max(round(m),300);%3*m
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
Q=X{1};
Qn=length(Q);
mineig_Q=eig(Q);
mineig_Q=min(mineig_Q);
Q=Q-mineig_Q*eye(Qn);
err=A*Q(:)-b;
err=norm(err,1);
% disp('err')
% disp(err)
% for i=1:6
%     if(err<0.5)
%         break
%     end
%     disp(i)
%     [obj,X,s,y,S,Z,ybar,v,info,runhist]= sdpnalplus(blk,At,C,b,[],[],[],[],[],opti,X,s,y,S,Z,ybar,v);
%     Q=X{1};
%     Qn=length(Q);
%     mineig_Q=eig(Q);
%     mineig_Q=min(mineig_Q);
%     Q=Q-mineig_Q*eye(Qn);
%     err=A*Q(:)-b;
%     err=norm(err,1);
%     disp('err')
%     disp(err)
%
% end
% [obj,X,s,y,S,Z,ybar,v,info,runhist]= sdpnalplus(blk,At,C,b,[],[],[],[],[],opti,X,s,y,S,Z,ybar,v);
lb=norm(A*Q(:)-b,1)+trace(Q)- length(Q)*lambda_min(Q);

end