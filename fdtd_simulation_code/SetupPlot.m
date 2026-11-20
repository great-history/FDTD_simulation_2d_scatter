function [I1] = SetupPlot(EzMatrix, yPosList, xPosList, lambda, monLeftIdx, monRightIdx, monUpperIdx, monLowerIdx)
    figure
    set(gcf, 'unit', 'inch', 'position', [3, 3, 7, 6],'doublebuffer', 'on') % figure
    set(gca, 'Position', [0.2, 0.12, 0.70, 0.75]) % axess
    I1 = imagesc(xPosList, yPosList, gather(EzMatrix), [-1.5,1.5]);
    % colormap('jet')
    colormap(parula(256));
    colorbar;
    hold on;
    viscircles([0 0], 150*10^(-9), 'Color', 'r', 'LineWidth', 1);
    hold on;
    xCorners = [xPosList(monLeftIdx), xPosList(monRightIdx), xPosList(monRightIdx), xPosList(monLeftIdx)];  % x坐标，比如按顺时针排列
    yCorners = [yPosList(monUpperIdx), yPosList(monUpperIdx), yPosList(monLowerIdx), yPosList(monLowerIdx)];  % y坐标
    plot([xCorners, xCorners(1)], [yCorners, yCorners(1)], 'r--', 'LineWidth', 1);
    
    title(['2D map of Ez ( 0 \% finished and $\lambda = $', num2str(lambda), ' nm)'], 'Interpreter','latex');
    xlabel('x-[m]', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex');
    ylabel('y-[m]', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex');

    drawnow;
end


%