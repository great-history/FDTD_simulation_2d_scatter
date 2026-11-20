%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% set up : source is the Modulated Gaussian Pulse (MPG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lambda, omega, numFreq, omegaList, ampEz, tau, t0, dy, ...
            yNum, yPosList, dx, xNum, xPosList, circLinIdx, dt, timeSteps, ...
                pmlWidth, sfwidth, tfXWidth, tfYWidth, monitorWidth, epsEnv, muEnv, epsMatrix, muMatrix, ...
                    sigmaEXProfile, sigmaEYProfile, sigmaMXProfile, sigmaMYProfile] = Dielectric_Cylinder_Setup_MGP()
    %% parameters
    mu0 = 4 * pi * 10^(-7); % vaccum mu0, unit : T·m/A
    eps0 = 8.8541878188 * 10^(-12); % vaccum eps0, unit : F·m^{-1}
    cLight = 1 / sqrt(mu0 * eps0); % c_light, unit : m/s
    nsi = 3.48; % refraction index of Silicon, unit : ??
    
    %% source : plane wave
    ampEz = 1;
    lambda = 600.0; % unit : nm   450 / 500 / 550 / 600 / 650 / 700 / 750
    tPeriod = (lambda * 10^(-9)) / cLight;
    omega = 2 * pi / tPeriod;
    omegaHi = 2 * pi / (450 * 10^(-9)) * cLight;
    omegaLo = 2 * pi / (750 * 10^(-9)) * cLight;
    tau = (2 * pi) / (omegaHi - omegaLo);
    t0 = 2.5 * tau;
    
    %% spatial resolution
    % the size of the Silicon Dielectric
    radius = 150 * 10^(-9); % unit : nm

    % dx = min((lambda / (nsi * 15))* 10^(-9), radius/10);
    % dx = (10 * 10^(-9)) * floor(dx * 10^(9) / 10); % make dx to be integer multiple of 10 nm
    % xLen = max(radius*2*10, 8 * lambda * 10^(-9));
    % xNum = ceil((xLen / dx) / 50) * 50 + 1; % number of cells along x direction, must be odd
    
    dx = min((lambda/nsi)/20*10^(-9), radius/30);
    xNum = 401; % 
    xPosList = (0:xNum - 1) * dx;
    xPosList = xPosList - mean(xPosList);
    % xLen = xNum * dx; % unit : nm, physical size along x direction
    
    dy = dx; % step size along y direction, unit : m
    yNum = xNum; % yNum = 1001; % number of cells along y direction, must be odd
    yPosList = (0:yNum - 1) * dy;
    yPosList = yPosList - mean(yPosList);
    % yLen = yNum * dy; % unit : nm, physical size along y direction

    % xCTRIdx 

    % create a circle at the center : Silicon Scatter
    [xMesh, yMesh] = meshgrid(xPosList, yPosList);
    circMask = ((xMesh).^2 + (yMesh).^2) <= (radius)^2;
    circLinIdx = find(circMask);
    
    %% time resolution
    sFactor = 1 / sqrt(2);
    % dt = min((sFactor * 0.9) * dx / (cLight * nsi), tPeriod / 100); % unit : s
    dt = min((sFactor * 0.9) * dx / cLight, tPeriod / 100); % unit : s
    timeSteps = 50 * round(tPeriod / dt);
    
    % frequency
    kminIdx = floor(omegaLo / (2 * pi) * (timeSteps * dt));
    kmaxIdx = ceil(omegaHi / (2 * pi) * (timeSteps * dt));
    numFreq = kmaxIdx - kminIdx + 1;
    omegaList = (2 * pi / timeSteps) * linspace(kminIdx, kmaxIdx, numFreq);

    %% grid
    pmlWidth = 25; % 75 the width of the PML layers on each side
    sfwidth = 100; % 200 the width of the scattered field region on each side
    tfXWidth = xNum - 2 * (pmlWidth + sfwidth);
    tfYWidth = yNum - 2 * (pmlWidth + sfwidth);
    % monitor ( on a large square enclosed the total field region )
    monitorWidth = xNum - 2 * pmlWidth - sfwidth; % the width of a monitor

    %% for a silicon (n = 3.48) dielectric cylinder, mu = mu0, & eps = eps0 * n^2
    epsEnv = eps0; muEnv = mu0; 
    epsScat = eps0 * nsi^2; muScat = mu0; 
    % epsScat = eps0; muScat = mu0;
    muPML = mu0; epsPML = eps0; % permittivity & permeability

    %% PML parameters
    Rf0 = 1e-8; mOrder = 6; % reflectance coefficient at incident angle = 0 / use a polynomial profile for conductivities
    
    %% get all medium properties
    [epsMatrix, muMatrix, sigmaEXProfile, sigmaEYProfile, sigmaMXProfile, sigmaMYProfile] = ...
                get_all_medium_properties(yNum, xNum, circLinIdx, epsEnv, muEnv, epsScat, muScat, epsPML, muPML, Rf0, mOrder, pmlWidth, dx);

end