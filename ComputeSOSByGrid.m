function [Q_min,lb_min]=ComputeSOSByGrid(Q,A,b,Maxiter)
for i=1:1
    Q_old=Q;
    g=Grid(Q,A,b);g=g+g';g=g/2;
    t=1e-6;
    beta=0.1;
    lb= norm(A*Q(:)-b,1)+trace(Q)- length(Q)*lambda_min(Q);
    Q0=Q;
    f0=norm(A*Q0(:)-b,1)+trace(Q0)- length(Q)*lambda_min(Q0);
    while(f0>lb-0.1*t*norm(g(:))^2)
        t=beta*t;
        Q0=Q+t*g;
        f0=norm(A*Q0(:)-b,1)+trace(Q0)- length(Q)*lambda_min(Q0);
    end
    D_Q=Q0-Q;
    Q=Q0;
    lb_min=f0;
    Q_min=Q;
end

%% looping
small_step_numbers=0;
for min_t=1e-16./(1:5)
    for i=1:Maxiter
        g_old=g;
        Q_old=Q;
        g=Grid(Q,A,b);g=g+g';g=g/2;
        %     disp(norm(t));
        D_g=-1*(g-g_old);
        t1= (D_Q(:)'*D_g(:))/(D_g(:)'*D_g(:));
        t2=(D_Q(:)'*D_Q(:))/(D_Q(:)'*D_g(:));
        t=max(t1,t2);
        t=t2;
        t=max(min_t,min(t,1/length(Q)));
        Q=Q+t*g;
        lb=norm(A*Q(:)-b,1)+trace(Q)- length(Q)*lambda_min(Q);
%         disp(vpa([t,lb,[norm(D_Q(:))^2,norm(D_g(:))^2],norm(g(:))],8))
        if lb<lb_min
            Q_min=Q;
            lb_min=norm(A*Q(:)-b,1)+trace(Q)- length(Q)*lambda_min(Q);
        end
        if norm(g(:))<1e-6
            lb_min
            return
        end
        if  norm(D_Q(:))>1e-6
            small_step_numbers=0;
        else
            small_step_numbers=small_step_numbers+1;
            if small_step_numbers>=10 ;       lb_min;     return;end
        end
        D_Q=Q-Q_old;
    end
end
lb_min
end


function [alpha,f0,Q0]=linesearch(Q,G,A,b,alpha,rho,c)
lb=norm(A*Q(:)-b,1)+trace(Q)- length(Q)*lambda_min(Q);
Q0=Q+alpha*G;
f0=norm(A*Q0(:)-b,1)+trace(Q0)- length(Q)*lambda_min(Q0);%f0=f(Q0)

%% W method:
while(f0> lb-c*alpha*norm(G(:))^2)
    alpha=alpha*rho;
    Q0=Q+alpha*G;
    f0=norm(A*Q0(:)-b,1)+trace(Q0)- length(Q)*lambda_min(Q0);
end
% disp(vpa([f0,f0-lb,norm(Q-Q0),norm(G)],5))
%% Backtrick
% beta=0.8;
% while(f0> lb-alpha/2*norm(G)^2)
%     alpha=alpha*beta;
%     Q0=Q+alpha*G;
%     f0=norm(A*Q0(:)-b,1)+trace(Q0)- length(Q)*lambda_min(Q0);
% end
%   disp(vpa([f0,f0-lb,norm(Q-Q0)],3))
%%
end

function [g,A,b]=Grid2(Q,A,b)
%min tr(Q)+||AQ(:)-b)||_1-len(Q)*min_eig(Q)
[a,~]=eig(Q);
n=length(Q);
u=a(:,1);
G1=-n*u*u';
y=A*Q(:)-b;
y2=A'*sign(y);
G2=reshape(y2,[n,n]);
G3=eye(n);G3(1)=0;
g=G1+G2+G3;
g=-g;
% g=-g;

end