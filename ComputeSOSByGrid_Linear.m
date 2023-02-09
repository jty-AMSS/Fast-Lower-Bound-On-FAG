function [Q_min,lb_min]=ComputeSOSByGrid_Linear(X0,A,b,Suppf_lb,Maxinneriter,tol)
if nargin <=4
    Maxinneriter=500;
end
if nargin <=5
    tol=1e-10;
end
Q=X0;
% g=Grid(Q,A,b);g=g+g';g=g/2;
% t=0.001;
% Q=Q+t*g;
lb_min=norm(A*Q(:)-b,1)+trace(Q)- length(Q)*lambda_min(Q);
%% looping
small_step_numbers=0;
for i=1:Maxinneriter
    g=Grid(Q,A,b);g=g+g';g=g/2;
    [alpha,lb_min,Q]=linesearch(Q,g,A,b,1,0.8,1/1000);
    if norm(g(:))<tol||alpha*norm(g(:))<1e-16
        break
    end
    if  alpha*norm(g(:))>1e-6
        small_step_numbers=0;
    else
        small_step_numbers=small_step_numbers+1;
        if small_step_numbers>=5;       lb_min;     break;end
    end
end
lb_min=lb_min-Suppf_lb;
Q_min=Q;
end



function [alpha,f0,Q0]=linesearch(Q,G,A,b,alpha,rho,c)
lb=norm(A*Q(:)-b,1)+trace(Q)- length(Q)*lambda_min(Q);
Q0=Q+alpha*G;
f0=norm(A*Q0(:)-b,1)+trace(Q0)- length(Q)*lambda_min(Q0);%f0=f(Q0)

%% W method:
% clc
while(f0> lb-c*alpha*norm(G(:))^2)
    alpha=alpha*rho;
    Q0=Q+alpha*G;
    f0=norm(A*Q0(:)-b,1)+trace(Q0)- length(Q)*lambda_min(Q0);
%     disp(vpa([f0-lb,alpha],4))
end
disp(vpa([f0,f0-lb,norm(Q-Q0),norm(G(:))],5))
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