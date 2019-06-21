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
%       x0 - 4 X 1 state initialization
%       F  - 2 X nt system input (force)
%       C  - 4 X 4 observation matrix
%       y_true - sampled observation
% output: 
%       pe - prediction error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pe = pe1_func(params,tspan,x0,C,F,y_true)
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
    [dim_obs,~] = size(C);
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
    state = zeros(4,length(tspan));
    ysol = zeros(dim_obs,length(tspan));
    state(:,1) = x0;
    for i = 1:1:length(tspan)-1
        state(:,i+1) = A*state(:,i)+B*[0;0;F(1,i);F(2,i)];
        ysol(:,i+1) = C*state(:,i+1);
    end    
    pe = norm(ysol-y_true,2);
end