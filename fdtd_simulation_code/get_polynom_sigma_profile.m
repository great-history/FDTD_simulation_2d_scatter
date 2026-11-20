%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% function to get polynomial profile of electric conductivity
%%%%% assume dx = dy (unit : m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigmaEXProfile, sigmaEYProfile] = get_polynom_sigma_profile(Rf0, mOrder, pmlWidth, refracIdx, dx)
    boundList = (1:1:pmlWidth);
    
    sigmaEXProfile = - ((mOrder + 1) * log(Rf0)) / (2 * refracIdx * pmlWidth * dx) * (boundList / pmlWidth).^mOrder;
    sigmaEYProfile = sigmaEXProfile; % they are the same in this case
end