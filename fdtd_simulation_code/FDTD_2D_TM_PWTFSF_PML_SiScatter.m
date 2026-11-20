%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% 2D FDTD TM mode ( Ez, Hx, Hy ) with plane wave source in TF/SF framework 
%%%%%% does not works
%%%%%% to do list : symmetry / only update tIdx / GPU acceleration / TE mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
%% parameters
% vaccum
mu0 = 4 * pi * 10^(-7); % unit : T·m/A
eps0 = 8.8541878188 * 10^(-12); % unit : F·m^{-1}
cLight = 1 / sqrt(mu0 * eps0); % unit : m/s

%% source : plane wave
% source : 1. ramped sinusoid / 2. unit-step sinusoid
% assume wave propagation direction is along y direction, Ez : sin(ky - \omega t), Hx = \sqrt{eps/mu} sin(ky - \omega t)
lambda = 750.0; % unit : nm   450 / 500 / 550 / 600 / 650 / 700 / 750
tPeriod = (lambda * 10^(-9)) / cLight;
omega = 2 * pi / tPeriod;
ampEz = 1; % unit : N/C
ampHx = sqrt(eps0 / mu0) * ampEz;
ki = omega / cLight; % unit : 1 / m

%% spatial discretization of the 2D computational domain
dx = 2.5e-9; % step size along x direction, unit : m
xNum = 1001; % number of cells along x direction
dy = dx; % step size along y direction, unit : m
yNum = 1001; % number of cells along y direction
radius = 150 * 10^(-9); % unit : nm
[yPosList, xPosList, circLinIdx] = create_grid_with_circle(dx, dy, xNum, yNum, radius);

%% time resolution
sFactor = 1 / sqrt(2);
dt = min(sFactor * dx / cLight * 0.99, tPeriod / 100); % unit : s
timeSteps = 15 * round(tPeriod / dt);

%% Index
pmlWidth = 75; % the width of the PML layers on each side
sfwidth = 200; % the width of the scattered field region on each side
tfXWidth = xNum - pmlWidth - sfwidth;
tfYWidth = yNum - pmlWidth - sfwidth;

% monitor ( on a large square enclosed the total field region )
monitorWidth = xNum - 2 * pmlWidth - sfwidth; % the width of a monitor

[ywdIdx, xwdIdx, pmlUpperIdx, pmlLowerIdx, pmlLeftIdx, pmlRightIdx, ycIdx, xcIdx, ...
    tsUpperIdx, tsLowerIdx, tsLeftIdx, tsRightIdx, monUpperIdx, monLowerIdx, monLeftIdx, monRightIdx] = ...
            get_indexes(yNum, xNum, pmlWidth, sfwidth, monitorWidth);
monitorWidth = monRightIdx - monLeftIdx + 1;
numMonitor = 4 * (monitorWidth - 1);

%% initialization
% for a silicon (n = 3.48) dielectric cylinder, mu = mu0, & eps = eps0 * n^2
epsEnv = eps0; muEnv = mu0; epsScat = eps0 * 3.48^2; muScat = mu0; muPML = mu0; epsPML = eps0; % permittivity & permeability
Rf0 = 1e-8; mOrder = 6; % reflectance coefficient at incident angle = 0 / use a polynomial profile for conductivities
[epsMatrix, muMatrix, sigmaEXProfile, sigmaEYProfile, sigmaMXProfile, sigmaMYProfile] = ...
                get_all_medium_properties(yNum, xNum, circLinIdx, epsEnv, muEnv, epsScat, muScat, epsPML, muPML, Rf0, mOrder, pmlWidth, dx);
% for 1D update
[IncHxEz, IncHyEz, IncEzHx, IncEzHy] = get_1D_update_coefficients(dt, dy, dx, epsEnv, muEnv);

% for 2D update coefficients for the whole region
CHxHx = ones(yNum, xNum);
CHxEz = (dt / dy) ./ muMatrix;

CHyHy = ones(yNum, xNum);
CHyEz = (dt / dx) ./ muMatrix;

CEzEz = ones(yNum, xNum);
CEzHx = (dt / dy) ./ epsMatrix;
CEzHy = (dt / dx) ./ epsMatrix;

[KEzEz, KEzHx, KEzHy, PEzxEzx, PEzxHy, KHyHy, KHyEz, PEzyEzy, PEzyHx, KHxHx, KHxEz] = ...
        get_2D_update_coefficients(yNum, xNum, pmlWidth, dt, dy, dx, epsMatrix, muMatrix, ...
                                sigmaEXProfile, sigmaMXProfile, sigmaEYProfile, sigmaMYProfile);
%% field components
[EzMatrix, HxMatrix, HyMatrix, EzIncList, HxIncList, HyIncList, EzxPML, EzyPML, ...
            tfsfUpperList, tfsfLowerList, tfsfLeftList, tfsfRightList, EzMonList, HxMonList, HyMonList] = ...
                intialize_field_components(yNum, xNum, timeSteps, pmlWidth, numMonitor, tfXWidth, tfYWidth);

% total scattering cross section
totSCSList = zeros(timeSteps+1, 1);

%% update loop
tic
figure
set(gcf, 'unit', 'inch', 'position', [3, 3, 13.5, 5],'doublebuffer', 'on') % figure
set(gca, 'Position', [0.2, 0.12, 0.70, 0.75]) % axess

for tIdx = 0:timeSteps
    %% turn on HxHy monitor (at time (n-1/2)*dt)
    HxMonList(1:monitorWidth) = HxMatrix(monUpperIdx, monLeftIdx:monRightIdx);
    HxMonList(monitorWidth+1:2*monitorWidth) = HxMatrix(monLowerIdx, monLeftIdx:monRightIdx);
    
    HyMonList(1:monitorWidth) = HyMatrix(monUpperIdx:monLowerIdx, monLeftIdx);
    HyMonList(monitorWidth+1:2*monitorWidth) = HyMatrix(monUpperIdx:monLowerIdx, monRightIdx);

    %% 1D update : assume n-3/2 & n-1 has obtained, source is located at tsUpperIndex - 1
    tnow = tIdx * dt;
    % Hx at time (n - 1/2) * dt : incident field
    HxIncList(1:end-1) = HxIncList(1:end-1) - IncHxEz * (EzIncList(2:end) - EzIncList(1:end-1));
    
    % Ez at time n * dt : incident field
    EzIncList(2:end) = EzIncList(2:end) - IncEzHx * (HxIncList(2:end) - HxIncList(1:end-1));
    % EzIncList(tsUpperIdx-1) = ampEz * sin(omega * tnow);
    EzIncList(1) = ampEz * sin(omega * tnow);
    
    %% 2D update : assume (n-1) & (n-1/2) has obtained
    %% here present a modified version of updating, which may take more time but can get exactly zero incident waves in the scattered region
    %% Ez at time n * dt : total field & scattered field ( need Hx,Hy at time (n-1/2) * dt )
    % update the total field region without TF/SF interface included
    % extract the SF/TF interface
    
    EzMatrix(ycIdx, xcIdx) = KEzEz(ycIdx, xcIdx) .* EzMatrix(ycIdx, xcIdx) ...
                            + KEzHy(ycIdx, xcIdx) .* (HyMatrix(ycIdx, xcIdx) - HyMatrix(ycIdx, xcIdx-1)) ...
                            - KEzHx(ycIdx, xcIdx) .* (HxMatrix(ycIdx, xcIdx) - HxMatrix(ycIdx-1, xcIdx));
    % special treatment at SF/TF interface ( need Hx & Hy at time (n - 1/2) * dt )
    
    % within the front face
    EzMatrix(tsUpperIdx, tsLeftIdx:tsRightIdx) = EzMatrix(tsUpperIdx, tsLeftIdx:tsRightIdx) + IncEzHx * HxIncList(1);

    % within the back face
    EzMatrix(tsLowerIdx, tsLeftIdx:tsRightIdx) = EzMatrix(tsLowerIdx, tsLeftIdx:tsRightIdx) - IncEzHx * HxIncList(tsLowerIdx - tsUpperIdx + 2);
    
    % within Left Face
    EzMatrix(tsUpperIdx:tsLowerIdx, tsLeftIdx) = EzMatrix(tsUpperIdx:tsLowerIdx, tsLeftIdx) ...
        - CEzHy(tsUpperIdx:tsLowerIdx, tsLeftIdx) .* HyIncList(2:tsLowerIdx-tsUpperIdx+2);

    % within right face
    EzMatrix(tsUpperIdx:tsLowerIdx, tsRightIdx) = EzMatrix(tsUpperIdx:tsLowerIdx, tsRightIdx) ...
        + CEzHy(tsUpperIdx:tsLowerIdx, tsRightIdx) .* HyIncList(2:tsLowerIdx-tsUpperIdx+2);
    
    %% Ez at time n * dt : PML layers
    % b. update Ezx, Ezy, Ez at PML layers
    % left PML layer
    EzxPML(:, 2:pmlWidth) = PEzxEzx(:, pmlWidth-1:-1:1) .* EzxPML(:, 2:pmlWidth) ...
                + PEzxHy(:, pmlWidth-1:-1:1) .* (HyMatrix(:, 2:pmlWidth) - HyMatrix(:, 1:pmlWidth-1));

    EzyPML(ycIdx, 1:pmlWidth) = EzyPML(ycIdx, 1:pmlWidth)...
                - CEzHx(1,1) * (HxMatrix(ycIdx, 1:pmlWidth) - HxMatrix(ycIdx-1, 1:pmlWidth));
    
    % right PML layer
    EzxPML(:, pmlRightIdx) = PEzxEzx(:, 1:pmlWidth) .* EzxPML(:, pmlRightIdx) ...
                + PEzxHy(:, 1:pmlWidth) .* (HyMatrix(:, pmlRightIdx) - HyMatrix(:, pmlRightIdx-1));

    EzyPML(ycIdx, pmlRightIdx) = EzyPML(ycIdx, pmlRightIdx) ...
                - CEzHx(1,1) * (HxMatrix(ycIdx, pmlRightIdx) - HxMatrix(ycIdx-1, pmlRightIdx));

    % upper PML layer
    EzxPML(pmlUpperIdx, xcIdx) = EzxPML(pmlUpperIdx, xcIdx)...
                + CEzHy(1,1) * (HyMatrix(pmlUpperIdx, xcIdx) - HyMatrix(pmlUpperIdx, xcIdx-1));

    EzyPML(2:pmlWidth, :) = PEzyEzy(pmlWidth-1:-1:1, :) .* EzyPML(2:pmlWidth, :) ...
                - PEzyHx(pmlWidth-1:-1:1, :) .* (HxMatrix(2:pmlWidth, :) - HxMatrix(1:pmlWidth-1, :));

    % lower PML layer
    EzxPML(pmlLowerIdx, xcIdx) = EzxPML(pmlLowerIdx, xcIdx) ...
                + CEzHy(1,1) * (HyMatrix(pmlLowerIdx, xcIdx) - HyMatrix(pmlLowerIdx, xcIdx-1));
    EzyPML(pmlLowerIdx, :) = PEzyEzy(1:pmlWidth, :) .* EzyPML(pmlLowerIdx, :) ...
                - PEzyHx(1:pmlWidth, :) .* (HxMatrix(pmlLowerIdx, :) - HxMatrix(pmlLowerIdx-1, :));

    % update Ez in the PML layer
    EzMatrix(pmlUpperIdx, :) = EzxPML(pmlUpperIdx, :) + EzyPML(pmlUpperIdx, :);
    EzMatrix(pmlLowerIdx, :) = EzxPML(pmlLowerIdx, :) + EzyPML(pmlLowerIdx, :);
    EzMatrix(ycIdx, pmlLeftIdx) = EzxPML(ycIdx, pmlLeftIdx) + EzyPML(ycIdx, pmlLeftIdx);
    EzMatrix(ycIdx, pmlRightIdx) = EzxPML(ycIdx, pmlRightIdx) + EzyPML(ycIdx, pmlRightIdx);
    
    %% Hx, Hy at time (n + 1/2) * dt : total field & scattered field ( need Ez at time n * dt )
    HxMatrix(ywdIdx, xwdIdx) = KHxHx(ywdIdx, xwdIdx) .* HxMatrix(ywdIdx, xwdIdx) ...
                            - KHxEz(ywdIdx, xwdIdx) .* (EzMatrix(ywdIdx+1, xwdIdx) - EzMatrix(ywdIdx, xwdIdx));
    
    HyMatrix(ywdIdx, xwdIdx) = KHyHy(ywdIdx, xwdIdx) .* HyMatrix(ywdIdx, xwdIdx) ...
                            + KHyEz(ywdIdx, xwdIdx) .* (EzMatrix(ywdIdx, xwdIdx+1) - EzMatrix(ywdIdx, xwdIdx));

    % special treatment at SF/TF interface ( need Ez at time n * dt )

    % outside the front face
    HxMatrix(tsUpperIdx-1, tsLeftIdx:tsRightIdx) = HxMatrix(tsUpperIdx-1, tsLeftIdx:tsRightIdx)...
        + CHxEz(tsUpperIdx-1, tsLeftIdx:tsRightIdx) * EzIncList(2);
    
    % outside the back face
    HxMatrix(tsLowerIdx, tsLeftIdx:tsRightIdx) = HxMatrix(tsLowerIdx, tsLeftIdx:tsRightIdx)...
        - CHxEz(tsLowerIdx, tsLeftIdx:tsRightIdx) .* EzIncList(tsLowerIdx - tsUpperIdx + 2);
    
    % outside the left face
    HyMatrix(tsUpperIdx:tsLowerIdx, tsLeftIdx-1) = HyMatrix(tsUpperIdx:tsLowerIdx, tsLeftIdx-1)...
        - CHyEz(tsUpperIdx:tsLowerIdx, tsLeftIdx-1) .* EzIncList(2:tsLowerIdx-tsUpperIdx+2);

    % outside the right face
    HyMatrix(tsUpperIdx:tsLowerIdx, tsRightIdx) = HyMatrix(tsUpperIdx:tsLowerIdx, tsRightIdx)...
        + CHyEz(tsUpperIdx:tsLowerIdx, tsRightIdx) .* EzIncList(2:tsLowerIdx-tsUpperIdx+2);
    
    %% turn on Ez monitor
    EzMonList(1:monitorWidth-1) = EzMatrix(monUpperIdx, monLeftIdx:monRightIdx-1);
    EzMonList(monitorWidth:2*monitorWidth-2) = EzMatrix(monUpperIdx:monLowerIdx-1, monRightIdx);
    EzMonList(2*monitorWidth-1:3*monitorWidth-3) = EzMatrix(monLowerIdx, monRightIdx:-1:monLeftIdx+1);
    EzMonList(3*monitorWidth-2:4*monitorWidth-4) = EzMatrix(monLowerIdx:-1:monUpperIdx+1, monLeftIdx);
    
    %% turn on HxHy monitor again
    HxMonList(1:monitorWidth) = (HxMonList(1:monitorWidth) + transpose(HxMatrix(monUpperIdx, monLeftIdx:monRightIdx))) / 2;
    HxMonList(monitorWidth+1:2*monitorWidth) = (HxMonList(monitorWidth+1:2*monitorWidth) + transpose(HxMatrix(monLowerIdx, monLeftIdx:monRightIdx))) / 2;
    
    HyMonList(1:monitorWidth) = (HyMonList(1:monitorWidth) + HyMatrix(monLowerIdx:-1:monUpperIdx, monLeftIdx)) / 2;
    HyMonList(monitorWidth+1:2*monitorWidth) = (HyMonList(monitorWidth+1:2*monitorWidth) + HyMatrix(monLowerIdx:-1:monUpperIdx, monRightIdx)) / 2;
    
    %% calculate total scattering cross section
    % corners
    totSCSList(tIdx + 1) = 1 / sqrt(2) * EzMonList(1) * (- HxMonList(1) - HyMonList(1)); % upper left corner
    totSCSList(tIdx + 1) = totSCSList(tIdx + 1) + 1 / sqrt(2) * EzMonList(monitorWidth) * (HxMonList(monitorWidth) - HyMonList(monitorWidth+1)); % upper right corner
    totSCSList(tIdx + 1) = totSCSList(tIdx + 1) + 1 / sqrt(2) * EzMonList(2*monitorWidth-1) * (HxMonList(2*monitorWidth) + HyMonList(2*monitorWidth)); % lower right corner
    totSCSList(tIdx + 1) = totSCSList(tIdx + 1) + 1 / sqrt(2) * EzMonList(2*monitorWidth-2) * (- HxMonList(monitorWidth+1) + HyMonList(monitorWidth)); % lower left corner
    % upper side and lower side
    totSCSList(tIdx + 1) = totSCSList(tIdx + 1) - dot(EzMonList(2:monitorWidth-1), HxMonList(2:monitorWidth-1)); % the upper side (normal vector : -y)
    totSCSList(tIdx + 1) = totSCSList(tIdx + 1) + dot(EzMonList(3*monitorWidth-3:-1:2*monitorWidth), HxMonList(monitorWidth+2:2*monitorWidth-1)); % the lower side (normal vector : +y)
    % left side and right side
    totSCSList(tIdx + 1) = totSCSList(tIdx + 1) + dot(EzMonList(4*monitorWidth-4:-1:3*monitorWidth-1), HyMonList(2:monitorWidth-1)); % the left side (normal vector : -x)
    totSCSList(tIdx + 1) = totSCSList(tIdx + 1) - dot(EzMonList(monitorWidth+1:2*monitorWidth-2), HyMonList(monitorWidth+2:2*monitorWidth-1)); % the right side (normal vector : +x)
    %% movie show
    if mod(tIdx, 20) == 1
        subplot(1,2,1)
        imagesc(xPosList, yPosList, EzMatrix, [-1.5,1.5]);
        colorbar;
        hold on;
        viscircles([0 0], 150*10^(-9), 'Color', 'r', 'LineWidth', 1);
        hold on;
        xCorners = [xPosList(monLeftIdx), xPosList(monRightIdx), xPosList(monRightIdx), xPosList(monLeftIdx)];  % x坐标，比如按顺时针排列
        yCorners = [yPosList(monUpperIdx), yPosList(monUpperIdx), yPosList(monLowerIdx), yPosList(monLowerIdx)];  % y坐标
        plot([xCorners, xCorners(1)], [yCorners, yCorners(1)], 'r--', 'LineWidth', 1);

        title(['2D map of Ez at time = ', num2str(tnow * 1e+15, '%.2f'), ' fs']);
        xlabel('x(in m)', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex');
        ylabel('y(in m)', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex')
    
        % subplot(1,2,2)
        % plot(EzIncList);
        % axis([0 timeSteps/5 -1.25, 1.25])
        % title(['1D wave', num2str(tnow * 1e+15, '%.2f'), ' fs']);
        % xlabel('y(in nm)', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex');
        % ylabel('Ez', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex')
        
        subplot(1,2,2)
        plot(EzMonList);
        axis([1 numMonitor -1.25, 1.25])
        title(['Ez on the monitor', num2str(tnow * 1e+15, '%.2f'), ' fs']);
        xlabel('x/dx', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex');
        ylabel('Ez', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex')

        getframe;
    end
end
toc

if ~exist('currentDateTime', 'var')
    currentDateTime = datestr(now, 'yyyymmdd_HHMMSS');
end
%% save img
fdtdPlot = figure();
imagesc(xPosList, yPosList, EzMatrix, [-1.5,1.5]);
colorbar;
hold on;
viscircles([0 0], 150*10^(-9), 'Color', 'r', 'LineWidth', 1);
hold on;
xCorners = [xPosList(monLeftIdx), xPosList(monRightIdx), xPosList(monRightIdx), xPosList(monLeftIdx)];  % x坐标，比如按顺时针排列
yCorners = [yPosList(monUpperIdx), yPosList(monUpperIdx), yPosList(monLowerIdx), yPosList(monLowerIdx)];  % y坐标
plot([xCorners, xCorners(1)], [yCorners, yCorners(1)], 'r--', 'LineWidth', 1);

title(['2D map of Ez at time = ', num2str(tnow * 1e+15, '%.2f'), ' ( ', num2str(lambda), ' nm) ']);
xlabel('x(in m)', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex');
ylabel('y(in m)', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Arial', 'Interpreter','latex')

folderPath = '.\fdtd_saved_image_folder';
fileName = [currentDateTime '_lambda_' num2str(lambda) 'nm_fdtd.pdf'];
% 构建完整的文件路径
fullFilePath = fullfile(folderPath, fileName);
exportgraphics(gcf, fullFilePath, 'ContentType', 'vector', 'BackgroundColor','none');
close(fdtdPlot)

%% save data
folderPath = '.\fdtd_saved_data_folder';
dataName = [currentDateTime '_lambda_' num2str(lambda) 'nm.mat'];
fullFilePath = fullfile(folderPath, dataName);
save(fullFilePath);