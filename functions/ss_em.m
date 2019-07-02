%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EM-algorithm for a state space model
% Jingxiao Liu
% June 2019
% Input: 
%       A    -- transition matrix
%       B  -- input matrix
%       C  -- observation matrix
%       Q  -- covariance of process noise
%       R  -- covariance of observation noise
%       Y  -- observation
%       x0 -- state at the first time stamp
% output: 
%       XT -- X[k|0:T]
%       PT -- P[k|0:T]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,C,Q,R,llh] = ss_em(A,C,Q,R,Y,x0,p0)
num = size(Y,2);
tol = 1e-6;
maxit = 1000;
llh = -inf(1,maxit);
for iter = 2:maxit
%     E-step
    [XT, PT, Exx, Ex1x, llh(iter)] = kalman_smoother(A,C,Q,R,Y,x0,p0);
    if abs(llh(iter)-llh(iter-1)) < tol*abs(llh(iter-1)) %check convergence
        break;
    end   
%     M-step 
    EXX = sum(Exx,3);
    EX1X = sum(Ex1x,3);
    A = EX1X/(EXX-Exx(:,:,num));
    Q = (EXX-Exx(:,:,1)-EX1X*A')/(num-1);
    Q = (Q+Q')/2;
%     C = Y*XT'/EXX;
    R = (Y*Y'-Y*XT'*C')/num;
    R = (R+R')/2;
    x0 = XT(:,1);
    p0 = PT(:,:,1);    
end
llh = llh(2:iter);
end
