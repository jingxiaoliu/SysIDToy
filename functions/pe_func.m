%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prediction error calculation: 2DoF ocsillator
% data
% Jingxiao Liu
% 5/30/2019
% |\/\/\/\O\/\/\/\O
%  k1 c1 m1 k2 c2 m2
% Input: 
%       params - system parameters
%       tspan  - time span
%       F  - 2 X nt system input (force)
%       y_true - sampled observation
% output: 
%       pe - prediction error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pe = pe_func(params,tspan,F,y_true)
    m1 = params(1);
    m2 = params(2);
    k1 = params(3);
    k2 = params(4);
    c1 = params(5);
    c2 = params(6);
    M=[m1 0;0 m2]; %mass matrix
    K=[k1+k2 -k2; -k2 k2]; %stiffness matrix
    Damp=[c1+c2 -c2; -c2 c2]; % damping matrix
    dt = tspan(2)-tspan(1); % delta time
    [dim_obs,~] = size(y_true);
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
    yp = zeros(dim_obs,length(tspan));
    for i = 1:1:length(tspan)-1
        yp(:,i+1) = A*y_true(:,i)+B*[0;0;F(1,i);F(2,i)];
    end    
    pe = norm(yp-y_true,2);
end