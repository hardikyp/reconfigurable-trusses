clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath("SequentialKinematics/");
addpath("TopologyOptimization/");
addpath("TraditionalTrusses/");
addpath("ModifiedTopOptTruss/");
addpath("Videos/");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run("WarrenTruss.m");

%% === IDENTIFY TRIANGULAR ELEMENTS =======================================
[processedData] = FindTriangles(inputData);

%% === LOCATION FOR NEW NODES ENABLING FLAT FOLDABILITY ===================
[flatFoldableTruss] = GenerateFFTruss(processedData);

%% === TARGET ANGLES AND SEARCH GRID ======================================
targetAnglesDeg = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 15, 20];
steps = 10;
% Search over pctAct in [0, 1]. Increase pctActRange(2) if a target angle
% is reported as unreachable for any unit.
pctActRange = [0, 1];
pctActSamples = unique(sort([ ...
    pctActRange(1), ...
    logspace(-8, log10(pctActRange(2)), 250), ...
    linspace(pctActRange(1), pctActRange(2), 250)]));

%% === SAMPLE pctAct -> ANGLE MAP =========================================
numSamples = numel(pctActSamples);
angleSamplesDeg = zeros(numSamples, 3);

for i = 1:numSamples
    angleSamplesDeg(i, :) = computeUnitAngles(flatFoldableTruss, ...
                                              pctActSamples(i), ...
                                              steps);
end

%% === FIND pctAct FOR EACH TARGET ANGLE ==================================
pctActResults = NaN(numel(targetAnglesDeg), 3);
for unitIdx = 1:3
    for targetIdx = 1:numel(targetAnglesDeg)
        pctActResults(targetIdx, unitIdx) = findPctActForTarget( ...
            flatFoldableTruss, ...
            pctActSamples, ...
            angleSamplesDeg(:, unitIdx), ...
            targetAnglesDeg(targetIdx), ...
            unitIdx, ...
            steps);
    end
end

resultsTable = table(targetAnglesDeg(:), ...
                     pctActResults(:, 1), ...
                     pctActResults(:, 2), ...
                     pctActResults(:, 3), ...
                     'VariableNames', ...
                     {'targetAngle_deg', ...
                      'pctAct_unit_1', ...
                      'pctAct_unit_2', ...
                      'pctAct_unit_3'});

reachableRangeTable = table((1:3)', ...
                            min(angleSamplesDeg, [], 1)', ...
                            max(angleSamplesDeg, [], 1)', ...
                            'VariableNames', ...
                            {'unit', ...
                             'minAngle_deg', ...
                             'maxAngle_deg'});

disp(resultsTable);
disp(reachableRangeTable);


function anglesDeg = computeUnitAngles(flatFoldableTruss, pctAct, steps)
    numDOF = length(flatFoldableTruss.newNodeNum);
    percentActuationDOF = ones(1, numDOF) * pctAct;
    kinematicsData = SequentialKinematics(flatFoldableTruss, steps, ...
                                          percentActuationDOF);

    coords = kinematicsData.kinematicsNodeLoc(:, :, end);
    unit_1_angle = 180 - abs(rad2deg(atan2(coords(2, 2) - coords(3, 2), ...
                                           coords(2, 1) - coords(3, 1)) - ...
                                     atan2(coords(4, 2) - coords(3, 2), ...
                                           coords(4, 1) - coords(3, 1))));
    unit_2_angle = 180 - abs(rad2deg(atan2(coords(4, 2) - coords(5, 2), ...
                                           coords(4, 1) - coords(5, 1)) - ...
                                     atan2(coords(7, 2) - coords(5, 2), ...
                                           coords(7, 1) - coords(5, 1))));
    unit_3_angle = 180 - abs(rad2deg(atan2(coords(7, 2) - coords(8, 2), ...
                                           coords(7, 1) - coords(8, 1)) - ...
                                     atan2(coords(9, 2) - coords(8, 2), ...
                                           coords(9, 1) - coords(8, 1))));
    anglesDeg = [unit_1_angle, unit_2_angle, unit_3_angle];
end


function pctAct = findPctActForTarget(flatFoldableTruss, pctActSamples, ...
                                      angleSamplesDeg, targetAngleDeg, ...
                                      unitIdx, steps)
    pctAct = NaN;

    minAngle = min(angleSamplesDeg);
    maxAngle = max(angleSamplesDeg);
    if targetAngleDeg < minAngle || targetAngleDeg > maxAngle
        return;
    end

    diffSamples = angleSamplesDeg - targetAngleDeg;
    exactIdx = find(abs(diffSamples) < 1e-12, 1);
    if ~isempty(exactIdx)
        pctAct = pctActSamples(exactIdx);
        return;
    end

    % Use a monotonic interpolation as the first estimate, then refine with
    % fzero if a sign-changing bracket exists in the sampled response.
    [sortedAngles, sortIdx] = sort(angleSamplesDeg);
    sortedPctAct = pctActSamples(sortIdx);
    [uniqueAngles, uniqueIdx] = unique(sortedAngles, 'stable');
    uniquePctAct = sortedPctAct(uniqueIdx);
    pctGuess = interp1(uniqueAngles, uniquePctAct, targetAngleDeg, ...
                       'pchip', NaN);

    bracketIdx = find(diffSamples(1:end-1) .* diffSamples(2:end) < 0);
    if isempty(bracketIdx)
        pctAct = pctGuess;
        return;
    end

    if numel(bracketIdx) > 1 && ~isnan(pctGuess)
        bracketCenters = 0.5 * (pctActSamples(bracketIdx) + ...
                                pctActSamples(bracketIdx + 1));
        [~, bestIdx] = min(abs(bracketCenters - pctGuess));
        bracketIdx = bracketIdx(bestIdx);
    else
        bracketIdx = bracketIdx(1);
    end

    objective = @(pct) computeSingleUnitAngle(flatFoldableTruss, pct, ...
                                              steps, unitIdx) - ...
                       targetAngleDeg;
    lowerBound = pctActSamples(bracketIdx);
    upperBound = pctActSamples(bracketIdx + 1);

    try
        pctAct = fzero(objective, [lowerBound, upperBound]);
    catch
        pctAct = pctGuess;
    end
end


function angleDeg = computeSingleUnitAngle(flatFoldableTruss, pctAct, ...
                                           steps, unitIdx)
    anglesDeg = computeUnitAngles(flatFoldableTruss, pctAct, steps);
    angleDeg = anglesDeg(unitIdx);
end
