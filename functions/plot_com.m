%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for plotting comparison
% Generator
% Jingxiao Liu
% June 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_com(t,y,ycom,ylabels,com_legend1,com_legend2)
num_signals = size(y,1);
for i = 1:num_signals
    subplot(num_signals,1,i)    
    plot(t,y(i,:),'b'); hold on
    plot(t,ycom(i,:),'g'); hold off
    xlabel('time');
    ylabel(ylabels{i});
    legend(com_legend1,com_legend2);
end