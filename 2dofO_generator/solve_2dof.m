%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for solving differential equaiton
% velocity
% Generator
% Jingxiao Liu
% 5/30/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt=solve_2dof(t,y,M,C,K,Ft,F)
f = [interp1(Ft,F(1,:),t);interp1(Ft,F(2,:),t)];
dydt = zeros(size(y));
dydt(1:2)= y(3:4);
dydt2 = inv(M)*(f.')' - (inv(M)*K*y(1:2))-(inv(M)*C*y(3:4));
dydt(3:4) = dydt2;