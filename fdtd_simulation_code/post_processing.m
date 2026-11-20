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

for kIdx = 1:numFreq
    EzMonList = EzOmegaList(:, kIdx);
    HxMonList = HxOmegaList(:, kIdx);
    HxMonList = conj(HxMonList);
    HyMonList = HyOmegaList(:, kIdx);
    HyMonList = conj(HyMonList);
    
    % upper left corner
    totSCS = 1 / sqrt(2) * EzMonList(1) * (- HxMonList(1) - HyMonList(1));
    % upper right corner
    totSCS = totSCS + 1 / sqrt(2) * EzMonList(monitorWidth) * (HxMonList(monitorWidth) - HyMonList(monitorWidth+1));
    % lower right corner
    totSCS = totSCS + 1 / sqrt(2) * EzMonList(2*monitorWidth-1) * (HxMonList(2*monitorWidth) + HyMonList(2*monitorWidth));
    % lower left corner
    totSCS = totSCS + 1 / sqrt(2) * EzMonList(2*monitorWidth-2) * (- HxMonList(monitorWidth+1) + HyMonList(monitorWidth));

    % the upper side (normal vector : -y)
    totSCS = totSCS - dot(EzMonList(2:monitorWidth-1), HxMonList(2:monitorWidth-1));
    % the lower side (normal vector : +y)
    totSCS = totSCS + dot(EzMonList(3*monitorWidth-3:-1:2*monitorWidth), HxMonList(monitorWidth+2:2*monitorWidth-1));
    % the left side (normal vector : -x)
    totSCS = totSCS + dot(EzMonList(4*monitorWidth-4:-1:3*monitorWidth-1), HyMonList(2:monitorWidth-1));
    % the right side (normal vector : +x)
    totSCS = totSCS - dot(EzMonList(monitorWidth+1:2*monitorWidth-2), HyMonList(monitorWidth+2:2*monitorWidth-1));
    
    % get the magnitude of the fourier component of source
    sourceFT = tau * sqrt(pi) / 2 * (exp(-((omegaList(kIdx)/dt - omega) * tau / 2)^2) + exp(-((omegaList(kIdx)/dt + omega) * tau / 2)^2));
    incIntensity = epsEnv * sourceFT^2 / 2 * cLight; % unit : W / m^2
    totSCSList(kIdx) = totSCS / incIntensity;
end

% plot
totSCSPlot = figure();
plot(omegaList, real(totSCSList * dx / incIntensity * 10^(18))); % integral element : dx, multiply by 10^(18) is to convert the unit from m^2 to nm^2
% axis([omegaList(1), omegaList(end), 0, 5*10^(11)])
title(['$\sigma_{tot}(t)$ ', ' ($\lambda = $', num2str(lambda), ' nm)'], 'Interpreter','latex');
xlabel('$\frac{t}{dt}$', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex');
ylabel('$\sigma_{tot}(t) - [nm^2]$', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex')

% folderPath = '.\fdtd_saved_image_folder';
% fileName = [currentDateTime '_lambda_' num2str(lambda) 'nm_totSCS.pdf'];
% % 构建完整的文件路径
% fullFilePath = fullfile(folderPath, fileName);
% exportgraphics(gcf, fullFilePath, 'ContentType', 'vector', 'BackgroundColor','none');
% close(totSCSPlot)

%% save data
folderPath = '.\fdtd_saved_data_folder';
dataName = [currentDateTime '_lambda_' num2str(lambda) 'nm.mat'];
fullFilePath = fullfile(folderPath, dataName);
save(fullFilePath, 'currentDateTime', 'lambda', 'dt', 'dx', 'dy', 'timeSteps', 'yNum', 'xNum', 'totSCSList', 'EzMatrix', 'HxMatrix', 'HyMatrix');