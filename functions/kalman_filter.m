%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kalman filter
% Jingxiao Liu
% June 2019
% Input: 
%       A  -- transition matrix
%       C  -- observation matrix
%       Q  -- covariance of process noise
%       R  -- covariance of observation noise
%       Y  -- observation
%       x0 -- state at the first time stamp
% output: 
%       Xk1k -- X[k+1|k]
%       Xkk  -- X[k|k]
%       Pk1k -- P[k+1|k]
%       Pkk  -- P[k|k]
%       llh  -- log-likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xk1k,Xkk,Pk1k,Pkk,llh] = kalman_filter(A,C,Q,R,Y,x0,p0)
xkk = x0;
pkk = p0;
num = size(Y,2);
Xk1k = zeros(4, num);
Xkk = zeros(4, num);
Pkk = zeros(4, 4, num);
Pk1k = zeros(4, 4, num);
logpdf = zeros(1,num);
for i = 1:size(Y,2)
    xk1k = A * xkk;
    pk1k = A * pkk * A' + Q;
    K = pk1k*C'*inv(R + C*pk1k*C');
    xkk = xk1k + K * (Y(:,i)-C*xk1k);
    pkk = pk1k - K*C*pk1k;
    Xk1k(:,i) = xk1k;
    Xkk(:,i) = xkk;
    Pk1k(:,:,i) = pk1k;
    Pkk(:,:,i) = pkk;
    
    U = chol(R + C*pk1k*C');
    logpdf(i) = - (0.5 * size(Y,1) * log (2 * pi)...
                + 0.5 * sum (((U') \ (Y(:,i) - C*xk1k)) .^ 2)...
                + sum (log (diag (U))));
end
llh = sum(logpdf);
end

