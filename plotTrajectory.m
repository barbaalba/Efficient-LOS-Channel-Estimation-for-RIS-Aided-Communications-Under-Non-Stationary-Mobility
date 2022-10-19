function plotTrajectory(x_t,y_t,Az,conf)
K = size(x_t,1);
T = size(x_t,2);
for k = 1:K
    figure('defaultAxesFontSize',20,'DefaultLineLineWidth', 2);
    fig = gcf;
    set(fig,'position',[60 50 1600 800]); % [left bottom width height]
    %tiledlayout(2,1);
    subplot(1,2,1);set(gca, 'ydir', 'reverse');xlabel('y');ylabel('x');
    grid on; xlim([-25,25]); ylim([-25,25]); hold on;
    subplot(1,2,2);xlabel('time');ylabel('azimuth');grid on;ylim([-90,90]);
    xlim([1,T]);hold on;
    for t = 2:T
    %nexttile(1);
    subplot(1,2,1);
    if strcmp(conf,'discrete')
        hold off;
        plot(y_t(k,t-1:t),x_t(k,t-1:t),'Color','b'); set(gca, 'ydir', 'reverse');xlabel('y');ylabel('x');xlim([-25,25]); ylim([-25,25]);grid on;
    else
        plot(y_t(k,1:t),x_t(k,1:t),'Color','b');
    end
    %nexttile(2);
    subplot(1,2,2);
    plot(rad2deg(Az(k,1:t)),'r');
    pause(0.1);
    end
end