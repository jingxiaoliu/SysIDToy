%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian LTI system: 2DoF ocsillator data 
% Generator
% Jingxiao Liu
% 5/30/2019
% |\/\/\/\O\/\/\/\O
%  k1 c1 m1 k2 c2 m2
% Input: 
%       m1 - Mass of the first DoF
%       m2 - Mass of the second DoF
%       k1 - Stiffness of the first DoF
%       k2 - Stiffness of the second DoF
%       c1 - Damping coeff of the first DoF
%       c2 - Damping coeff of the second DoF
%       T  - Time duration of the output signal
%       nt - Number of time stamp
%       x0 - 4 X 1 state initialization
%       F  - 2 X nt system input (force)
%       C  - 4 X 4 observation matrix
%       promu - 4 X 1 mean vector of process 
%               noise
%       prosi - 4 X 4 covariance of process 
%               noise
%       obsmu - 4 X 1 mean vector of observation 
%               noise
%       obssi - 4 X 4 covariance of obserbation 
%               noise ~~~size change
% output: 
%       ysol1 - Observation [x1,x2,v1,v2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ysol = glti_2dofO_generator(m1,m2,k1,k2,c1,c2,T,nt,x0,F,C,...
                                        promu,prosi,obsmu,obssi,rs)
    rng(rs) %random seed
    M=[m1 0;0 m2]; %mass matrix
    K=[k1+k2 -k2; -k2 k2]; %stiffness matrix
    Damp=[c1+c2 -c2; -c2 c2]; % damping matrix
    tspan=linspace(0,T,nt);
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
        pro_error = mvnrnd(promu,prosi);
        obs_error = mvnrnd(obsmu,obssi);
        state(:,i+1) = A*state(:,i)+B*[0;0;F(1,i);F(2,i)]+pro_error';
        ysol(:,i+1) = C*state(:,i+1)+obs_error';
    end    
end