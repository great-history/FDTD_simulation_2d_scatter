if ~exist('currentDateTime', 'var')
    currentDateTime = datestr(now, 'yyyymmdd_HHMMSS');
end

%% save img of Ez
fdtdPlot = figure();
set(gcf, 'unit', 'inch', 'position', [2, 5, 15, 3.5],'doublebuffer', 'on') % figure
set(gca, 'Position', [0.1, 0.1, 0.80, 0.85]) % axess

colormap('jet')

subplot(1,3,1)
imagesc(xPosList, yPosList, EzMatrix, [-1.5,1.5]);
colorbar;
hold on;
viscircles([0 0], 150*10^(-9), 'Color', 'r', 'LineWidth', 1);
hold on;
xCorners = [xPosList(monLeftIdx), xPosList(monRightIdx), xPosList(monRightIdx), xPosList(monLeftIdx)];  % x坐标，比如按顺时针排列
yCorners = [yPosList(monUpperIdx), yPosList(monUpperIdx), yPosList(monLowerIdx), yPosList(monLowerIdx)];  % y坐标
plot([xCorners, xCorners(1)], [yCorners, yCorners(1)], 'r--', 'LineWidth', 1);

title(['2D map of Ez', ' ($\lambda = $', num2str(lambda), ' nm)'], 'Interpreter','latex');
xlabel('x(in m)', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex');
ylabel('y(in m)', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex')

subplot(1,3,2)
imagesc(xPosList, yPosList, HxMatrix, [-1.5/500,1.5/500]);
colorbar;
title(['2D map of Hx', ' ($\lambda = $', num2str(lambda), ' nm)'], 'Interpreter','latex');
xlabel('x(in m)', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex');
ylabel('y(in m)', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex')

subplot(1,3,3)
imagesc(xPosList, yPosList, HyMatrix, [-1.5/500,1.5/500]);
colorbar;
title(['2D map of Hy', ' ($\lambda = $', num2str(lambda), ' nm)'], 'Interpreter','latex');
xlabel('x(in m)', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex');
ylabel('y(in m)', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex')

folderPath = '.\fdtd_saved_image_folder';
fileName = [currentDateTime '_lambda_' num2str(lambda) 'nm_fdtd.pdf'];
% 构建完整的文件路径
fullFilePath = fullfile(folderPath, fileName);
exportgraphics(gcf, fullFilePath, 'ContentType', 'vector', 'BackgroundColor','none');
close(fdtdPlot)

%% save img of total scattering cross section
% get total scattering cross section
cLight = 1 / sqrt(muEnv * epsEnv);
Iinc = 1 / 2 * epsEnv * ampEz^2 * cLight;
crossSecList = totSCSList;
crossSecList = crossSecList * (dx) / Iinc; % unit : m^2
crossSecList = crossSecList * (10^(9))^2; % unit : nm^2
maxCrossSec = max(crossSecList);
maxCrossSec = ceil(maxCrossSec / 10^(12)) * 10^(12);
% plot
totSCSPlot = figure();
plot((1:timeSteps)*dt*10^(15), crossSecList); % integral element : dx, multiply by 10^(18) is to convert the unit from m^2 to nm^2
axis([0 timeSteps*dt*10^(15) 0 1.5*10^(12)])
title(['$\sigma_{tot}(t)$ ', ' ($\lambda = $', num2str(lambda), ' nm)'], 'Interpreter','latex');
xlabel('$t - [fs]$', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex');
ylabel('$\sigma_{tot}(t) - [nm^2]$', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex')

folderPath = '.\fdtd_saved_image_folder';
fileName = [currentDateTime '_lambda_' num2str(lambda) 'nm_totSCS.pdf'];
% 构建完整的文件路径
fullFilePath = fullfile(folderPath, fileName);
exportgraphics(gcf, fullFilePath, 'ContentType', 'vector', 'BackgroundColor','none');
close(totSCSPlot)

%% save data
folderPath = '.\fdtd_saved_data_folder';
dataName = [currentDateTime '_lambda_' num2str(lambda) 'nm.mat'];
fullFilePath = fullfile(folderPath, dataName);
save(fullFilePath, 'currentDateTime', 'lambda', 'dt', 'dx', 'dy', 'timeSteps', 'yNum', 'xNum', 'totSCSList', 'EzMatrix', 'HxMatrix', 'HyMatrix');