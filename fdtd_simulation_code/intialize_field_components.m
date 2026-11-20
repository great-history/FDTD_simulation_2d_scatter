function [EzMatrix, HxMatrix, HyMatrix, EzIncList, HxIncList, HyIncList, EzxPML, EzyPML, ...
            tfsfUpperList, tfsfLowerList, tfsfLeftList, tfsfRightList, EzMonList, HxMonList, HyMonList] = ...
                intialize_field_components(yNum, xNum, timeSteps, pmlWidth, monitorWidth, tfXWidth, tfYWidth)
    % for 2D update in whole region
    EzMatrix = zeros(yNum, xNum);
    HxMatrix = zeros(yNum, xNum);
    HyMatrix = zeros(yNum, xNum);
    
    % for special treatment at TF/SF interface, used for Ez & Hx & Hy
    tfsfUpperList = zeros(1, tfXWidth);
    tfsfLowerList = zeros(1, tfXWidth);
    tfsfLeftList = zeros(1, tfYWidth);
    tfsfRightList = zeros(1, tfYWidth);

    % for 1D update
    EzIncList = zeros(timeSteps+1, 1);
    HxIncList = zeros(timeSteps+1, 1);
    HyIncList = zeros(timeSteps+1, 1);
    
    % for 2D update in PML layers
    EzxPML = zeros(yNum, xNum);
    EzyPML = zeros(yNum, xNum);
    % EzxLeftPML = zeros(yNum, pmlWidth);
    % EzxRightPML = zeros(yNum, xNum);
    % EzyUpperPML = zeros(pmlWidth, xNum);
    % EzyLowerPML = zeros(pmlWidth, xNum);

    % for monitor
    numMonitor = 4 * (monitorWidth - 1);
    EzMonList = zeros(numMonitor, 1);
    HxMonList = zeros(2 * monitorWidth, 1);
    HyMonList = zeros(2 * monitorWidth, 1);
end

