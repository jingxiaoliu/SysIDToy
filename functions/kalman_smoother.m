%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kalman smoother
% Jingxiao Liu
% June 2019
% Input: 
%       A  -- transition matrix
%       B  -- input matrix
%       C  -- observation matrix
%       Q  -- covariance of process noise
%       R  -- covariance of observation noise
%       Y  -- observation
%       x0 -- state at the first time stamp
% output: 
%       XT   -- X[k|0:T]
%       PT   -- P[k|0:T]
%       Exx  -- E[x[k]x[k]^T]
%       Ex1x -- E[x[k+1]x[k]^T]
%       llh  -- log-likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [XT, PT, Exx, Ex1x, llh] = kalman_smoother(A,C,Q,R,Y,x0,p0)
[Xk1k,Xkk,Pk1k,Pkk,llh] = kalman_filter(A,C,Q,R,Y,x0,p0);
XT = zeros(4, size(Y,2));
PT = zeros(4, 4, size(Y,2));
XT(:,size(Y,2)) = Xkk(:,size(Y,2));
PT(:,:,size(Y,2)) = Pkk(:,:,size(Y,2));
Exx = zeros(4, 4, size(Y,2));
Ex1x = zeros(4, 4, size(Y,2)-1);
for k=(size(Xk1k,2)-1):-1:1
    L(:,:,k) = Pkk(:,:,k) * A' * inv(Pk1k(:,:,k+1)+diag(ones(1,4)*1e-7));
    XT(:,k) = Xkk(:,k) + L(:,:,k) * ...
        (XT(:,k+1) - Xk1k(:,k+1));
    PT(:,:,k) = Pkk(:,:,k) + L(:,:,k) * ...
        (PT(:,:,k+1) - Pk1k(:,:,k+1)) * L(:,:,k)';
    Exx(:,:,k) = PT(:,:,k) + XT(:,k)*XT(:,k)';
    Ex1x(:,:,k) = PT(:,:,k+1)*L(:,:,k')' + XT(:,k+1)*XT(:,k)';
end
end
