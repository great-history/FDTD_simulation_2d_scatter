%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% epsEnv : the permittivity of environment ( free space )
%%%%% muEnv : the permeability of environment ( free space )
%%%%% epsScat : the permittivity of scatter
%%%%% muScat : the permeability of scatter
%%%%% epsPML : the permittivity of PML layers
%%%%% muPML : the permeability of PML layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KEzEz, KEzHx, KEzHy, PEzxEzx, PEzxHy, KHyHy, KHyEz, PEzyEzy, PEzyHx, KHxHx, KHxEz] = ...
                get_2D_update_coefficients(yNum, xNum, pmlWidth, dt, dy, dx, epsMatrix, muMatrix, ...
                                            sigmaEXProfile, sigmaMXProfile, sigmaEYProfile, sigmaMYProfile)
    %% update coefficients without conductivity
    KEzEz = ones(yNum, xNum);
    KEzHx = (dt / dy) ./ epsMatrix;
    KEzHy = (dt / dx) ./ epsMatrix;

    %% then calculate update coefficients for 2D update ( PML / SF / TF )
    % coefficients for PML layers for wave propagating along X direction (Left & Right)
    sigmaEXMatrix = zeros(yNum, pmlWidth); % E:electric / X:wave propagate along x direction
    sigmaMXMatrix = zeros(yNum, pmlWidth); % M : magnetic / X : wave propagate along x direction
    for yIdx = 1:yNum
        sigmaEXMatrix(yIdx, :) = sigmaEXProfile; % from small to large
        sigmaMXMatrix(yIdx, :) = sigmaMXProfile;
    end
    
    % P coefficients only for PML layers (Left & Right), K coefficients for the whole domain
    epsPMLXMatrix = epsMatrix(:, 1:pmlWidth);
    PEzxEzx = (epsPMLXMatrix - 0.5 * dt * sigmaEXMatrix) ./ (epsPMLXMatrix + 0.5 * dt * sigmaEXMatrix); % X:wave propagate along x direction
    PEzxHy = (dt / dx) ./ (epsPMLXMatrix + 0.5 * dt * sigmaEXMatrix);
    
    muPMLXMatrix = muMatrix(:, 1:pmlWidth);
    PHyHy = (muPMLXMatrix - 0.5 * dt * sigmaMXMatrix) ./ (muPMLXMatrix + 0.5 * dt * sigmaMXMatrix);
    KHyHy = ones(yNum, xNum); % update coefficients for the whole domain
    KHyHy(:, pmlWidth:-1:1) = PHyHy(:, 1:pmlWidth);
    KHyHy(:, xNum - pmlWidth + 1:xNum) = PHyHy(:, 1:pmlWidth);
    
    PHyEz = (dt / dx) ./ (muPMLXMatrix + 0.5 * dt * sigmaMXMatrix);
    KHyEz = (dt / dx) ./ muMatrix; % update coefficients for the whole domain
    KHyEz(:, pmlWidth:-1:1) = PHyEz(:, 1:pmlWidth);
    KHyEz(:, xNum - pmlWidth + 1:xNum) = PHyEz(:, 1:pmlWidth);
    
    % P coefficients for PML layers for wave propagating along Y direction, K coefficients for the whole domain
    sigmaEYMatrix = zeros(pmlWidth, xNum); % E : electric / Y :wave propagate along y direction
    sigmaMYMatrix = zeros(pmlWidth, xNum); % M : magnetic / Y : wave propagate along y direction
    for xIdx = 1:xNum
        sigmaEYMatrix(:, xIdx) = sigmaEYProfile; % from small to large
        sigmaMYMatrix(:, xIdx) = sigmaMYProfile;
    end
    
    % coefficients only for PML layers (upper & lower)
    epsPMLYMatrix = epsMatrix(1:pmlWidth, :);
    muPMLYMatrix = muMatrix(1:pmlWidth, :);
    PEzyEzy = (epsPMLYMatrix - 0.5 * dt * sigmaEYMatrix) ./ (epsPMLYMatrix + 0.5 * dt * sigmaEYMatrix); % Y:wave propagate along y direction
    PEzyHx = (dt / (dy)) ./ (epsPMLYMatrix + 0.5 * dt * sigmaEYMatrix);
    
    PHxHx = (muPMLYMatrix - 0.5 * dt * sigmaMYMatrix) ./ (muPMLYMatrix + 0.5 * dt * sigmaMYMatrix); % use the sigmaMY rather than sigmaMX ! ! !
    PHxEz = (dt / (dy)) ./ (muPMLYMatrix + 0.5 * dt * sigmaMYMatrix);
    
    KHxHx = ones(yNum, xNum); % update coefficients for the whole domain
    KHxHx(pmlWidth:-1:1, :) = PHxHx(1:pmlWidth, :);
    KHxHx(yNum - pmlWidth + 1:yNum, :) = PHxHx(1:pmlWidth, :);
    
    KHxEz = (dt / dy) ./ muMatrix;
    KHxEz(pmlWidth:-1:1, :) = PHxEz(1:pmlWidth, :);
    KHxEz(yNum - pmlWidth + 1:yNum, :) = PHxEz(1:pmlWidth, :);
end