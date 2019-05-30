%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for plotting displacement and
% velocity
% Generator
% Jingxiao Liu
% 5/30/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_dv(t,y)
subplot(2,1,1)
plot(t,y(:,1),'b'); hold on
plot(t,y(:,2),'g'); hold off
xlabel('time (s)');
ylabel('Displacement (m)');
legend('x1','x2');
subplot(2,1,2)
plot(t,y(:,3),'r'); hold on
plot(t,y(:,4),'k'); hold off
xlabel('time (s)');
ylabel('Velocity (m/s)');
legend('v1','v2');