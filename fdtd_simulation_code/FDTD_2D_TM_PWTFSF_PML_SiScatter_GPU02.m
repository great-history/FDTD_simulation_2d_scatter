%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% 2D FDTD TM mode ( Ez, Hx, Hy ) with plane wave source in TF/SF framework 
%%%%%% does not works
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lambda, omega, ampEz, dy, yNum, yPosList, dx, xNum, xPosList, circLinIdx, dt, timeSteps, ...
            pmlWidth, sfwidth, tfXWidth, tfYWidth, monitorWidth, epsEnv, muEnv, epsMatrix, muMatrix, ...
                sigmaEXProfile, sigmaEYProfile, sigmaMXProfile, sigmaMYProfile] = Dielectric_Cylinder_Setup();

%% Index
[ywdIdx, xwdIdx, pmlUpperIdx, pmlLowerIdx, pmlLeftIdx, pmlRightIdx, ycIdx, xcIdx, ...
    tsUpperIdx, tsLowerIdx, tsLeftIdx, tsRightIdx, monUpperIdx, monLowerIdx, monLeftIdx, monRightIdx, minSFXSteps, minPMLSteps] = ...
        get_indexes(yNum, xNum, pmlWidth, sfwidth, monitorWidth, circLinIdx);
monitorWidth = monRightIdx - monLeftIdx + 1;
numMonitor = 4 * (monitorWidth - 1);

%% update coefficients
% for 1D update
[IncHxEz, IncHyEz, IncEzHx, IncEzHy] = get_1D_update_coefficients(dt, dy, dx, epsEnv, muEnv);

% for 2D update coefficients for the whole region
[KEzEz, KEzHx, KEzHy, PEzxEzx, PEzxHy, KHyHy, KHyEz, PEzyEzy, PEzyHx, KHxHx, KHxEz] = ...
        get_2D_update_coefficients(yNum, xNum, pmlWidth, dt, dy, dx, epsMatrix, muMatrix, ...
                                sigmaEXProfile, sigmaMXProfile, sigmaEYProfile, sigmaMYProfile);

%% field components
[EzMatrix, HxMatrix, HyMatrix, EzIncList, HxIncList, HyIncList, EzxPML, EzyPML, ...
            tfsfUpperList, tfsfLowerList, tfsfLeftList, tfsfRightList, EzMonList, HxMonList, HyMonList] = ...
                intialize_field_components(yNum, xNum, timeSteps, pmlWidth, numMonitor, tfXWidth, tfYWidth);

% total scattering cross section
totSCSList = zeros(timeSteps - yNum, 1);

%% update loop
vars2GPU();
tic
figure
set(gcf, 'unit', 'inch', 'position', [3, 3, 13.5, 5],'doublebuffer', 'on') % figure
set(gca, 'Position', [0.2, 0.12, 0.70, 0.75]) % axess

for tIdx = 1:timeSteps
    %% 1D update : assume n-3/2 & n-1 has obtained, source is located at tsUpperIndex - 1
    tnow = (tIdx - 1) * dt;
    if tIdx < timeSteps
        % Hx at time (n - 1/2) * dt : incident field
        HxIncList(1:tIdx) = HxIncList(1:tIdx) - IncHxEz * (EzIncList(2:tIdx+1) - EzIncList(1:tIdx));
        
        % Ez at time n * dt : incident field
        EzIncList(2:tIdx+1) = EzIncList(2:tIdx+1) - IncEzHx * (HxIncList(2:tIdx+1) - HxIncList(1:tIdx));
        % EzIncList(tsUpperIdx-1) = ampEz * sin(omega * tnow);
        EzIncList(1) = ampEz * sin(omega * tnow);
    else
        % Hx at time (n - 1/2) * dt : incident field
        HxIncList(1:end-1) = HxIncList(1:end-1) - IncHxEz * (EzIncList(2:end) - EzIncList(1:end-1));
        
        % Ez at time n * dt : incident field
        EzIncList(2:end) = EzIncList(2:end) - IncEzHx * (HxIncList(2:end) - HxIncList(1:end-1));
        % EzIncList(tsUpperIdx-1) = ampEz * sin(omega * tnow);
        EzIncList(1) = ampEz * sin(omega * tnow);
    end
    
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
        - KEzHy(tsUpperIdx:tsLowerIdx, tsLeftIdx) .* HyIncList(2:tsLowerIdx-tsUpperIdx+2);

    % within right face
    EzMatrix(tsUpperIdx:tsLowerIdx, tsRightIdx) = EzMatrix(tsUpperIdx:tsLowerIdx, tsRightIdx) ...
        + KEzHy(tsUpperIdx:tsLowerIdx, tsRightIdx) .* HyIncList(2:tsLowerIdx-tsUpperIdx+2);
    
    %% Ez at time n * dt : PML layers
    % b. update Ezx, Ezy, Ez at PML layers
    % left PML layer
    EzxPML(:, 2:pmlWidth) = PEzxEzx(:, pmlWidth-1:-1:1) .* EzxPML(:, 2:pmlWidth) ...
                + PEzxHy(:, pmlWidth-1:-1:1) .* (HyMatrix(:, 2:pmlWidth) - HyMatrix(:, 1:pmlWidth-1));

    EzyPML(ycIdx, 1:pmlWidth) = EzyPML(ycIdx, 1:pmlWidth)...
                - KEzHx(1,1) * (HxMatrix(ycIdx, 1:pmlWidth) - HxMatrix(ycIdx-1, 1:pmlWidth));
    
    % right PML layer
    EzxPML(:, pmlRightIdx) = PEzxEzx(:, 1:pmlWidth) .* EzxPML(:, pmlRightIdx) ...
                + PEzxHy(:, 1:pmlWidth) .* (HyMatrix(:, pmlRightIdx) - HyMatrix(:, pmlRightIdx-1));

    EzyPML(ycIdx, pmlRightIdx) = EzyPML(ycIdx, pmlRightIdx) ...
                - KEzHx(1,1) * (HxMatrix(ycIdx, pmlRightIdx) - HxMatrix(ycIdx-1, pmlRightIdx));

    % upper PML layer
    EzxPML(pmlUpperIdx, xcIdx) = EzxPML(pmlUpperIdx, xcIdx)...
                + KEzHy(1,1) * (HyMatrix(pmlUpperIdx, xcIdx) - HyMatrix(pmlUpperIdx, xcIdx-1));

    EzyPML(2:pmlWidth, :) = PEzyEzy(pmlWidth-1:-1:1, :) .* EzyPML(2:pmlWidth, :) ...
                - PEzyHx(pmlWidth-1:-1:1, :) .* (HxMatrix(2:pmlWidth, :) - HxMatrix(1:pmlWidth-1, :));

    % lower PML layer
    EzxPML(pmlLowerIdx, xcIdx) = EzxPML(pmlLowerIdx, xcIdx) ...
                + KEzHy(1,1) * (HyMatrix(pmlLowerIdx, xcIdx) - HyMatrix(pmlLowerIdx, xcIdx-1));
    EzyPML(pmlLowerIdx, :) = PEzyEzy(1:pmlWidth, :) .* EzyPML(pmlLowerIdx, :) ...
                - PEzyHx(1:pmlWidth, :) .* (HxMatrix(pmlLowerIdx, :) - HxMatrix(pmlLowerIdx-1, :));

    % update Ez in the PML layer
    EzMatrix(pmlUpperIdx, :) = EzxPML(pmlUpperIdx, :) + EzyPML(pmlUpperIdx, :);
    EzMatrix(pmlLowerIdx, :) = EzxPML(pmlLowerIdx, :) + EzyPML(pmlLowerIdx, :);
    EzMatrix(ycIdx, pmlLeftIdx) = EzxPML(ycIdx, pmlLeftIdx) + EzyPML(ycIdx, pmlLeftIdx);
    EzMatrix(ycIdx, pmlRightIdx) = EzxPML(ycIdx, pmlRightIdx) + EzyPML(ycIdx, pmlRightIdx);
    
    %% turn on HxHy monitor (at time (n-1/2)*dt) % with spatial average
    if tIdx > yNum
        HxMonList(1:monitorWidth) = (HxMatrix(monUpperIdx-1, monLeftIdx:monRightIdx) + HxMatrix(monUpperIdx, monLeftIdx:monRightIdx))/2;
        HxMonList(monitorWidth+1:2*monitorWidth) = (HxMatrix(monLowerIdx, monLeftIdx:monRightIdx) + HxMatrix(monLowerIdx-1, monLeftIdx:monRightIdx))/2;
        
        HyMonList(1:monitorWidth) = (HyMatrix(monUpperIdx:monLowerIdx, monLeftIdx) + HyMatrix(monUpperIdx:monLowerIdx, monLeftIdx-1))/2;
        HyMonList(monitorWidth+1:2*monitorWidth) = (HyMatrix(monUpperIdx:monLowerIdx, monRightIdx) + HyMatrix(monUpperIdx:monLowerIdx, monRightIdx-1))/2;
    end

    %% Hx, Hy at time (n + 1/2) * dt : total field & scattered field ( need Ez at time n * dt )
    HxMatrix(ywdIdx, xwdIdx) = KHxHx(ywdIdx, xwdIdx) .* HxMatrix(ywdIdx, xwdIdx) ...
                            - KHxEz(ywdIdx, xwdIdx) .* (EzMatrix(ywdIdx+1, xwdIdx) - EzMatrix(ywdIdx, xwdIdx));
    
    HyMatrix(ywdIdx, xwdIdx) = KHyHy(ywdIdx, xwdIdx) .* HyMatrix(ywdIdx, xwdIdx) ...
                            + KHyEz(ywdIdx, xwdIdx) .* (EzMatrix(ywdIdx, xwdIdx+1) - EzMatrix(ywdIdx, xwdIdx));
    
    % special treatment at SF/TF interface ( need Ez at time n * dt )

    % outside the front face
    HxMatrix(tsUpperIdx-1, tsLeftIdx:tsRightIdx) = HxMatrix(tsUpperIdx-1, tsLeftIdx:tsRightIdx)...
        + KHxEz(tsUpperIdx-1, tsLeftIdx:tsRightIdx) * EzIncList(2);
    
    % outside the back face
    HxMatrix(tsLowerIdx, tsLeftIdx:tsRightIdx) = HxMatrix(tsLowerIdx, tsLeftIdx:tsRightIdx)...
        - KHxEz(tsLowerIdx, tsLeftIdx:tsRightIdx) .* EzIncList(tsLowerIdx - tsUpperIdx + 2);
    
    % outside the left face
    HyMatrix(tsUpperIdx:tsLowerIdx, tsLeftIdx-1) = HyMatrix(tsUpperIdx:tsLowerIdx, tsLeftIdx-1)...
        - KHyEz(tsUpperIdx:tsLowerIdx, tsLeftIdx-1) .* EzIncList(2:tsLowerIdx-tsUpperIdx+2);

    % outside the right face
    HyMatrix(tsUpperIdx:tsLowerIdx, tsRightIdx) = HyMatrix(tsUpperIdx:tsLowerIdx, tsRightIdx)...
        + KHyEz(tsUpperIdx:tsLowerIdx, tsRightIdx) .* EzIncList(2:tsLowerIdx-tsUpperIdx+2);
    
    if tIdx > yNum
        %% turn on Ez monitor
        EzMonList(1:monitorWidth-1) = EzMatrix(monUpperIdx, monLeftIdx:monRightIdx-1);
        EzMonList(monitorWidth:2*monitorWidth-2) = EzMatrix(monUpperIdx:monLowerIdx-1, monRightIdx);
        EzMonList(2*monitorWidth-1:3*monitorWidth-3) = EzMatrix(monLowerIdx, monRightIdx:-1:monLeftIdx+1);
        EzMonList(3*monitorWidth-2:4*monitorWidth-4) = EzMatrix(monLowerIdx:-1:monUpperIdx+1, monLeftIdx);
        
        %% turn on HxHy monitor again % with time average & spatial average
        HxMonList(1:monitorWidth) = HxMonList(1:monitorWidth) + transpose((HxMatrix(monUpperIdx-1, monLeftIdx:monRightIdx) + HxMatrix(monUpperIdx, monLeftIdx:monRightIdx))/2);
        HxMonList(1:monitorWidth) = HxMonList(1:monitorWidth) / 2;
        HxMonList(monitorWidth+1:2*monitorWidth) = HxMonList(monitorWidth+1:2*monitorWidth) + transpose((HxMatrix(monLowerIdx, monLeftIdx:monRightIdx) + HxMatrix(monLowerIdx-1, monLeftIdx:monRightIdx))/2);
        HxMonList(monitorWidth+1:2*monitorWidth) = HxMonList(monitorWidth+1:2*monitorWidth) / 2;
        
        HyMonList(1:monitorWidth) = HyMonList(1:monitorWidth) + (HyMatrix(monUpperIdx:monLowerIdx, monLeftIdx) + HyMatrix(monUpperIdx:monLowerIdx, monLeftIdx-1))/2;
        HyMonList(1:monitorWidth) = HyMonList(1:monitorWidth) / 2;
        HyMonList(monitorWidth+1:2*monitorWidth) = HyMonList(monitorWidth+1:2*monitorWidth) + (HyMatrix(monUpperIdx:monLowerIdx, monRightIdx) + HyMatrix(monUpperIdx:monLowerIdx, monRightIdx-1))/2;
        HyMonList(monitorWidth+1:2*monitorWidth) = HyMonList(monitorWidth+1:2*monitorWidth) / 2;
    
        %% calculate total scattering cross section
        totSCSList(tIdx) = 1 / sqrt(2) * EzMonList(1) * (- HxMonList(1) - HyMonList(1)); % upper left corner
        totSCSList(tIdx) = totSCSList(tIdx) + 1 / sqrt(2) * EzMonList(monitorWidth) * (HxMonList(monitorWidth) - HyMonList(monitorWidth+1)); % upper right corner
        totSCSList(tIdx) = totSCSList(tIdx) + 1 / sqrt(2) * EzMonList(2*monitorWidth-1) * (HxMonList(2*monitorWidth) + HyMonList(2*monitorWidth)); % lower right corner
        totSCSList(tIdx) = totSCSList(tIdx) + 1 / sqrt(2) * EzMonList(2*monitorWidth-2) * (- HxMonList(monitorWidth+1) + HyMonList(monitorWidth)); % lower left corner
        
        totSCSList(tIdx) = totSCSList(tIdx) - dot(EzMonList(2:monitorWidth-1), HxMonList(2:monitorWidth-1)); % the upper side (normal vector : -y)
        totSCSList(tIdx) = totSCSList(tIdx) + dot(EzMonList(3*monitorWidth-3:-1:2*monitorWidth), HxMonList(monitorWidth+2:2*monitorWidth-1)); % the lower side (normal vector : +y)
        
        totSCSList(tIdx) = totSCSList(tIdx) + dot(EzMonList(4*monitorWidth-4:-1:3*monitorWidth-1), HyMonList(2:monitorWidth-1)); % the left side (normal vector : -x)
        totSCSList(tIdx) = totSCSList(tIdx) - dot(EzMonList(monitorWidth+1:2*monitorWidth-2), HyMonList(monitorWidth+2:2*monitorWidth-1)); % the right side (normal vector : +x)
    end

    %% movie show
    % if t == 0
    %     [I1] = SetupPlot(Ez,Ne,x,y,dx,dy);
    % elseif mod(t,1) == 0
    %     E = abs(Ez(:,:,round(Nz/2))); E = 10*log10(E); 
    %     E(isinf(E)|isnan(E)) = 10*log10(eps(1+Single));
    %     I1.CData = gather(E); 
    %     drawnow;
    % end
    
    if mod(tIdx, 30) == 1
        subplot(1,2,1)
        imagesc(xPosList, yPosList, EzMatrix, [-1.5,1.5]);
        colormap('jet')
        colorbar;
        hold on;
        viscircles([0 0], 150*10^(-9), 'Color', 'r', 'LineWidth', 1);
        hold on;
        xCorners = [xPosList(monLeftIdx), xPosList(monRightIdx), xPosList(monRightIdx), xPosList(monLeftIdx)];  % x坐标，比如按顺时针排列
        yCorners = [yPosList(monUpperIdx), yPosList(monUpperIdx), yPosList(monLowerIdx), yPosList(monLowerIdx)];  % y坐标
        plot([xCorners, xCorners(1)], [yCorners, yCorners(1)], 'r--', 'LineWidth', 1);

        title(['2D map of Ez at time = ', num2str(tnow * 1e+15, '%.2f'), ' fs ($\lambda = $', num2str(lambda), ' nm)'], 'Interpreter','latex');
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