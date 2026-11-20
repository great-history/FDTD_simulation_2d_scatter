function [epsMatrix, muMatrix, sigmaEXProfile, sigmaEYProfile, sigmaMXProfile, sigmaMYProfile] = ...
                get_all_medium_properties(yNum, xNum, circLinIdx, epsEnv, muEnv, epsScat, muScat, epsPML, muPML, Rf0, mOrder, pmlWidth, dx)
    %% first calculate the material properties
    % assume isotropic at each point
    epsMatrix = epsEnv * ones(yNum, xNum);
    epsMatrix(circLinIdx) = epsScat;
    
    muMatrix = muEnv * ones(yNum, xNum);
    muMatrix(circLinIdx) = muScat;
    
    %% calculation of electrical & magnetic conductivity
    % electric and magnetic conductivity are not isotropic at the boundary
    refracIdx = sqrt((epsPML * muPML)/(epsEnv * muEnv)); % index of refraction, unit : dimensionless, with respect to environment which can be not vaccum
    [sigmaEXProfile, sigmaEYProfile] = get_polynom_sigma_profile(Rf0, mOrder, pmlWidth, refracIdx, dx);
    sigmaMXProfile = sigmaEXProfile * muPML / epsPML; % M : magnetic / X : wave propagate along x direction
    sigmaMYProfile = sigmaEYProfile * muPML / epsPML; % M : magnetic / Y : wave propagate along y direction
end