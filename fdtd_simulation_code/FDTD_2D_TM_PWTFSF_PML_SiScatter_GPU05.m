%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% 2D FDTD TM mode ( Ez, Hx, Hy ) with plane wave source in TF/SF framework 
%%%%%% using symmetry to reduce the computational time
%%%%%% using different types of sources : Continuous Plane Wave / Modulated Gaussian Pulse / Differential Gaussian Pulse
%%%%%% using fourier transformation to get total scattering cross section as a function of frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sourceType : 1 : Continuous Plane Wave / 2 : Modulated Gaussian Pulse / 3 : Differential Gaussian Pulse
sourceType = 2;
switch sourceType
    case 1
        [lambda, omega, ampEz, ...
            dy, yNum, yPosList, dx, xNum, xPosList, circLinIdx, dt, timeSteps, ...
                pmlWidth, sfwidth, tfXWidth, tfYWidth, monitorWidth, epsEnv, muEnv, epsMatrix, muMatrix, ...
                    sigmaEXProfile, sigmaEYProfile, sigmaMXProfile, sigmaMYProfile] = Dielectric_Cylinder_Setup();
        
        % fourier transformation
        omegaList = [0, omega, 2*omega];
    case 2
        [lambda, omega, numFreq, omegaList, ampEz, tau, t0, ...
            dy, yNum, yPosList, dx, xNum, xPosList, circLinIdx, dt, timeSteps, ...
                pmlWidth, sfwidth, tfXWidth, tfYWidth, monitorWidth, epsEnv, muEnv, epsMatrix, muMatrix, ...
                    sigmaEXProfile, sigmaEYProfile, sigmaMXProfile, sigmaMYProfile] = Dielectric_Cylinder_Setup_MGP();
        
    case 3
        [lambda, omega, ampEz, tau, t0, ... 
            dy, yNum, yPosList, dx, xNum, xPosList, circLinIdx, dt, timeSteps, ...
                pmlWidth, sfwidth, tfXWidth, tfYWidth, monitorWidth, epsEnv, muEnv, epsMatrix, muMatrix, ...
                    sigmaEXProfile, sigmaEYProfile, sigmaMXProfile, sigmaMYProfile] = Dielectric_Cylinder_Setup_DGP();
        
end


%% Index
[ywdIdx, xwdIdx, pmlUpperIdx, pmlLowerIdx, pmlLeftIdx, pmlRightIdx, ycIdx, xcIdx, ...
    tsUpperIdx, tsLowerIdx, tsLeftIdx, tsRightIdx, monUpperIdx, monLowerIdx, monLeftIdx, monRightIdx, minSFXSteps, minPMLSteps] = ...
        get_indexes(yNum, xNum, pmlWidth, sfwidth, monitorWidth, circLinIdx);
monitorWidth = monRightIdx - monLeftIdx + 1;
numMonitor = 4 * (monitorWidth - 1);

% x index of center of the domain
oxIdx = round(xNum/2);
xcLeftHalfIdx = xcIdx(1):oxIdx;
xcRightHalfIdx = oxIdx+1:xcIdx(end);
xwdLeftHalfIdx = xwdIdx(1):oxIdx-1;
xwdRightHalfIdx = xwdLeftHalfIdx + (oxIdx - 1);

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
                intialize_field_components(yNum, xNum, timeSteps, pmlWidth, monitorWidth, tfXWidth, tfYWidth);

% total scattering cross section
phaseFac = 1;
totSCSList = zeros(numFreq, 1); % total scattering cross section as a function of frequency omega
EzOmegaList = zeros(numMonitor, numFreq);
HxOmegaList = zeros(2 * monitorWidth, numFreq);
HyOmegaList = zeros(2 * monitorWidth, numFreq);

%% update loop
vars2GPU();

% setup plot
[I1] = SetupPlot(EzMatrix, yPosList, xPosList, lambda, monLeftIdx, monRightIdx, monUpperIdx, monLowerIdx);

tic
for tIdx = 1:timeSteps
    %% 1D update : assume n-3/2 & n-1 has obtained, source is located at tsUpperIndex - 1
    % Hx at time (n - 1/2) * dt : incident field
    HxIncList(1:tIdx) = HxIncList(1:tIdx) - IncHxEz * (EzIncList(2:tIdx+1) - EzIncList(1:tIdx));
    
    % Ez at time n * dt : incident field
    EzIncList(2:tIdx+1) = EzIncList(2:tIdx+1) - IncEzHx * (HxIncList(2:tIdx+1) - HxIncList(1:tIdx));
    
    tnow = (tIdx - 1) * dt;
    switch sourceType
        case 1
            EzIncList(1) = ampEz * sin(omega * tnow);
        case 2
            EzIncList(1) = ampEz * cos(omega * (tnow - t0)) * exp(-(tnow - t0)^2 / tau^2);
        case 3
            EzIncList(1) = - ampEz * (tnow - t0) / tau * exp(-(tnow - t0)^2 / tau^2);
    end
    
    %% 2D update : assume (n-1) & (n-1/2) has obtained
    %% here present a modified version of updating, which may take more time but can get exactly zero incident waves in the scattered region
    %% Ez at time n * dt : total field & scattered field ( need Hx,Hy at time (n-1/2) * dt )
    % update the total field region without TF/SF interface included
    % extract the SF/TF interface
    
    EzMatrix(ycIdx, xcLeftHalfIdx) = KEzEz(ycIdx, xcLeftHalfIdx) .* EzMatrix(ycIdx, xcLeftHalfIdx) ...
                            + KEzHy(ycIdx, xcLeftHalfIdx) .* (HyMatrix(ycIdx, xcLeftHalfIdx) - HyMatrix(ycIdx, xcLeftHalfIdx-1)) ...
                            - KEzHx(ycIdx, xcLeftHalfIdx) .* (HxMatrix(ycIdx, xcLeftHalfIdx) - HxMatrix(ycIdx-1, xcLeftHalfIdx));
    
    % special treatment at SF/TF interface ( need Hx & Hy at time (n - 1/2) * dt )
    
    % within the front face
    EzMatrix(tsUpperIdx, tsLeftIdx:tsRightIdx) = EzMatrix(tsUpperIdx, tsLeftIdx:tsRightIdx) + IncEzHx * HxIncList(1);
    
    % within the back face
    EzMatrix(tsLowerIdx, tsLeftIdx:tsRightIdx) = EzMatrix(tsLowerIdx, tsLeftIdx:tsRightIdx) - IncEzHx * HxIncList(tsLowerIdx - tsUpperIdx + 2);
    
    % within Left Face
    EzMatrix(tsUpperIdx:tsLowerIdx, tsLeftIdx) = EzMatrix(tsUpperIdx:tsLowerIdx, tsLeftIdx) ...
        - KEzHy(tsUpperIdx:tsLowerIdx, tsLeftIdx) .* HyIncList(2:tsLowerIdx-tsUpperIdx+2);
    
    % within right face
    % EzMatrix(tsUpperIdx:tsLowerIdx, tsRightIdx) = EzMatrix(tsUpperIdx:tsLowerIdx, tsRightIdx) ...
    %     + KEzHy(tsUpperIdx:tsLowerIdx, tsRightIdx) .* HyIncList(2:tsLowerIdx-tsUpperIdx+2);
    % by symmetry
    % EzMatrix(tsUpperIdx:tsLowerIdx, tsRightIdx) = EzMatrix(tsUpperIdx:tsLowerIdx, tsLeftIdx);
    
    % by symmetry, we have
    EzMatrix(ycIdx, xcRightHalfIdx) = EzMatrix(ycIdx, xcLeftHalfIdx(end-1:-1:1));

    %% Ez at time n * dt : PML layers
    % b. update Ezx, Ezy, Ez at PML layers
    % left PML layer
    EzxPML(:, 2:pmlWidth) = PEzxEzx(:, pmlWidth-1:-1:1) .* EzxPML(:, 2:pmlWidth) ...
                + PEzxHy(:, pmlWidth-1:-1:1) .* (HyMatrix(:, 2:pmlWidth) - HyMatrix(:, 1:pmlWidth-1));

    EzyPML(ycIdx, pmlLeftIdx) = EzyPML(ycIdx, pmlLeftIdx)...
                - KEzHx(1,1) * (HxMatrix(ycIdx, pmlLeftIdx) - HxMatrix(ycIdx-1, pmlLeftIdx));
    
    % right PML layer
    % EzxPML(:, pmlRightIdx) = PEzxEzx(:, 1:pmlWidth) .* EzxPML(:, pmlRightIdx) ...
    %             + PEzxHy(:, 1:pmlWidth) .* (HyMatrix(:, pmlRightIdx) - HyMatrix(:, pmlRightIdx-1));
    % 
    % EzyPML(ycIdx, pmlRightIdx) = EzyPML(ycIdx, pmlRightIdx) ...
    %             - KEzHx(1,1) * (HxMatrix(ycIdx, pmlRightIdx) - HxMatrix(ycIdx-1, pmlRightIdx));
    % by symmetry
    EzxPML(:, pmlRightIdx) = EzxPML(:, pmlLeftIdx(end:-1:1));
    EzyPML(:, pmlRightIdx) = EzyPML(:, pmlLeftIdx(end:-1:1));
    
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
    % EzMatrix(ycIdx, pmlRightIdx) = EzxPML(ycIdx, pmlRightIdx) + EzyPML(ycIdx, pmlRightIdx);
    EzMatrix(ycIdx, pmlRightIdx) = EzMatrix(ycIdx, pmlLeftIdx(end:-1:1)); % by symmetry

    %% turn on HxHy monitor (at time (n-1/2)*dt) % with spatial average
    HxMonList(1:monitorWidth) = (HxMatrix(monUpperIdx-1, monLeftIdx:monRightIdx) + HxMatrix(monUpperIdx, monLeftIdx:monRightIdx))/2;
    HxMonList(monitorWidth+1:2*monitorWidth) = (HxMatrix(monLowerIdx, monLeftIdx:monRightIdx) + HxMatrix(monLowerIdx-1, monLeftIdx:monRightIdx))/2;

    HyMonList(1:monitorWidth) = (HyMatrix(monUpperIdx:monLowerIdx, monLeftIdx) + HyMatrix(monUpperIdx:monLowerIdx, monLeftIdx-1))/2;
    HyMonList(monitorWidth+1:2*monitorWidth) = (HyMatrix(monUpperIdx:monLowerIdx, monRightIdx) + HyMatrix(monUpperIdx:monLowerIdx, monRightIdx-1))/2;

    %% Hx, Hy at time (n + 1/2) * dt : total field & scattered field ( need Ez at time n * dt )
    HxMatrix(ywdIdx, xwdLeftHalfIdx) = KHxHx(ywdIdx, xwdLeftHalfIdx) .* HxMatrix(ywdIdx, xwdLeftHalfIdx) ...
                            - KHxEz(ywdIdx, xwdLeftHalfIdx) .* (EzMatrix(ywdIdx+1, xwdLeftHalfIdx) - EzMatrix(ywdIdx, xwdLeftHalfIdx));
    % by symmetry
    HxMatrix(ywdIdx, xwdRightHalfIdx) = HxMatrix(ywdIdx, xwdLeftHalfIdx(end:-1:1));

    HyMatrix(ywdIdx, xwdLeftHalfIdx) = KHyHy(ywdIdx, xwdLeftHalfIdx) .* HyMatrix(ywdIdx, xwdLeftHalfIdx) ...
                            + KHyEz(ywdIdx, xwdLeftHalfIdx) .* (EzMatrix(ywdIdx, xwdLeftHalfIdx+1) - EzMatrix(ywdIdx, xwdLeftHalfIdx));
    % by symmetry
    HyMatrix(ywdIdx, xwdRightHalfIdx) = - HyMatrix(ywdIdx, xwdLeftHalfIdx(end:-1:1));

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
    
    %% turn on Ez monitor (at time n*dt) 
    EzMonList(1:monitorWidth-1) = EzMatrix(monUpperIdx, monLeftIdx:monRightIdx-1);
    EzMonList(monitorWidth:2*monitorWidth-2) = EzMatrix(monUpperIdx:monLowerIdx-1, monRightIdx);
    EzMonList(2*monitorWidth-1:3*monitorWidth-3) = EzMatrix(monLowerIdx, monRightIdx:-1:monLeftIdx+1);
    EzMonList(3*monitorWidth-2:4*monitorWidth-4) = EzMatrix(monLowerIdx:-1:monUpperIdx+1, monLeftIdx);

    %% turn on HxHy monitor again  (at time (n+1/2)*dt)  % with time average & spatial average
    HxMonList(1:monitorWidth) = HxMonList(1:monitorWidth) + transpose((HxMatrix(monUpperIdx-1, monLeftIdx:monRightIdx) + HxMatrix(monUpperIdx, monLeftIdx:monRightIdx))/2);
    HxMonList(1:monitorWidth) = HxMonList(1:monitorWidth) / 2;
    HxMonList(monitorWidth+1:2*monitorWidth) = HxMonList(monitorWidth+1:2*monitorWidth) + transpose((HxMatrix(monLowerIdx, monLeftIdx:monRightIdx) + HxMatrix(monLowerIdx-1, monLeftIdx:monRightIdx))/2);
    HxMonList(monitorWidth+1:2*monitorWidth) = HxMonList(monitorWidth+1:2*monitorWidth) / 2;

    HyMonList(1:monitorWidth) = HyMonList(1:monitorWidth) + (HyMatrix(monUpperIdx:monLowerIdx, monLeftIdx) + HyMatrix(monUpperIdx:monLowerIdx, monLeftIdx-1))/2;
    HyMonList(1:monitorWidth) = HyMonList(1:monitorWidth) / 2;
    HyMonList(monitorWidth+1:2*monitorWidth) = HyMonList(monitorWidth+1:2*monitorWidth) + (HyMatrix(monUpperIdx:monLowerIdx, monRightIdx) + HyMatrix(monUpperIdx:monLowerIdx, monRightIdx-1))/2;
    HyMonList(monitorWidth+1:2*monitorWidth) = HyMonList(monitorWidth+1:2*monitorWidth) / 2;
    
    %% calculate total scattering cross section
    for kIdx = 1:numFreq
        EzOmegaList(:, kIdx) = EzOmegaList(:, kIdx) + EzMonList * exp(-1j * omegaList(kIdx) * tIdx);
        HxOmegaList(:, kIdx) = HxOmegaList(:, kIdx) + HxMonList * exp(-1j * omegaList(kIdx) * tIdx);
        HyOmegaList(:, kIdx) = HyOmegaList(:, kIdx) + HyMonList * exp(-1j * omegaList(kIdx) * tIdx);
    end
    
    %% movie show
    if mod(tIdx, 50) == 1
        I1.CData = gather(EzMatrix);
        title(['2D map of Ez ( ', num2str(tIdx/timeSteps * 1e+2, '%.2f'), ' \% finished and $\lambda = $', num2str(lambda), ' nm)'], 'Interpreter','latex');
        drawnow;
    end
end
toc