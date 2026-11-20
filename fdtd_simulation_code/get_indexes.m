% pmlwidth : the width of the PML layers on each side
% sfwidth : the width of the scattered field region on each side
function [ywdIdx, xwdIdx, pmlUpperIdx, pmlLowerIdx, pmlLeftIdx, pmlRightIdx, ycIdx, xcIdx, ...
            tsUpperIdx, tsLowerIdx, tsLeftIdx, tsRightIdx, monUpperIdx, monLowerIdx, monLeftIdx, monRightIdx, minSFXSteps, minPMLSteps] = ...
                                                        get_indexes(yNum, xNum, pmlWidth, sfwidth, monitorWidth, circLinIdx)
    % whole domain indexes
    xwdIdx = 1 : xNum - 1;
    ywdIdx = 1 : yNum - 1;
    
    % PML layer indexes
    pmlLeftIdx = 1:pmlWidth;
    pmlRightIdx = xNum - pmlWidth + 1:xNum;
    pmlUpperIdx = 1:pmlWidth;
    pmlLowerIdx = yNum - pmlWidth + 1:yNum;
    
    % TF/SF framework
    % central region : includes both scattered field region and total field region
    ycIdx = pmlWidth + 1 : yNum - pmlWidth;
    xcIdx = pmlWidth + 1 : xNum - pmlWidth;
    
    % the TF/SF interface index
    tsLeftIdx = pmlWidth + sfwidth + 1; % this point is within the total field region
    tsRightIdx = xNum - sfwidth - pmlWidth; % this point is within the total field region
    tsUpperIdx = pmlWidth + sfwidth + 1; % this point is within the total field region
    tsLowerIdx = yNum - sfwidth - pmlWidth; % this point is within the total field region
    
    oxIdx = round(xNum/2);
    oyIdx = round(yNum/2);
    monUpperIdx = oyIdx - ceil(monitorWidth/2);
    monLowerIdx = oyIdx + ceil(monitorWidth/2);
    monLeftIdx = oxIdx - ceil(monitorWidth/2);
    monRightIdx = oxIdx + ceil(monitorWidth/2);
    
    % calculate the minimal steps the wave reach the left or right TF/SF boundary & left or right SF/PML boundary
    circUpperIdx = 1;
    for yIdx = 1:yNum
        if any(circLinIdx(yIdx, :))
            circUpperIdx = yIdx;
            break
        end
    end

    minSFXSteps = (circUpperIdx - tsUpperIdx + 1) + floor((tsRightIdx - tsLeftIdx + 1)/2); % the minimal steps reach the left / right TF/SF boundary
    minPMLSteps = (circUpperIdx - tsUpperIdx + 1) + floor((tsRightIdx - tsLeftIdx + 1)/2);
    minPMLSteps = min(minPMLSteps, tsLowerIdx - tsUpperIdx + 1);

    % minPMLStepsList = zeros(1,4); % upper / lower / left / right
    % minSFStepsList = zeros(1,4); % upper / lower / left / right
end