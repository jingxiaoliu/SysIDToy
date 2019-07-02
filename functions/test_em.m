clc;
clear;
close all;
addpath functions/
T = 10; %length of time duration
nt = 500; %number of time stamps
Fs = nt/T;
F = csvread(['../data/ambient.csv']); %load data
X = csvread(['../data/2dof_P6.csv']); %load data
tspan=linspace(0,T,nt);
m1=5; %mass 1 [kg]
m2=2; %mass 2 [kg]
k1=200; %spring 1 [N/m]
k2=100; %spring 2 [N/m]
c1=10; % damping coeff 1
c2=10; % damping coeff 2

% system matrices
M=[m1 0;0 m2]; %mass matrix
K=[k1+k2 -k2; -k2 k2]; %stiffness matrix
Damp=[c1+c2 -c2; -c2 c2]; % damping matrix
dt = tspan(2)-tspan(1); % delta time
C = [1,0,0,0;0,0,1,0];
Ac = [0,0,1,0;
     0,0,0,1;
     -K(1,1)/M(1,1),-K(1,2)/M(1,1),-Damp(1,1)/M(1,1),-Damp(1,2)/M(1,1);
     -K(2,1)/M(2,2),-K(2,2)/M(2,2),-Damp(2,1)/M(2,2),-Damp(2,2)/M(2,2)];
A = Ac*dt;
A = expm(A);
Bc = [0,0,0,0;
    0,0,0,0;
    0,0,1/M(1,1),0;
    0,0,0,1/M(2,2)];
B = inv(Ac)*(A-eye(4))*Bc;    
Y = C*X;
Q = diag([1e-6,1e-6,1e-6,1e-6]); % initialize process noise covariance
R = diag([1e-4,1e-4]); % initialize observation noise covariance
x0 = X(:,1);
p0 = Q;

num = size(Y,2);
tol = 1e-4;
maxit = 1000;
llh = -inf(1,maxit);
for iter = 2:maxit
%     E-step
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
    llh(iter) = sum(logpdf);
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
    if abs(llh(iter)-llh(iter-1)) < tol*abs(llh(iter-1)) %check convergence
        break;
    end   
%     M-step 
    EXX = sum(Exx,3);
    EX1X = sum(Ex1x,3);
    A = EX1X/(EXX-Exx(:,:,num));
    Q = (EXX-Exx(:,:,1)-EX1X*A')/(num-1);
    Q = (Q+Q')/2;
    C = Y*XT'/EXX;
    R = (Y*Y'-Y*XT'*C')/num;
    R = (R+R')/2;
    x0 = XT(:,1);
    p0 = PT(:,:,1);    
end
llh = llh(2:iter);