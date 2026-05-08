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

%% === KINEMATIC ANALYSIS =================================================
numDOF = length(flatFoldableTruss.newNodeNum);
percentActuationDOF = [0.9 , 0.86 , 0.9];
percentActuationDOF(1) = 7.8565e-05 ;%0.076607;
percentActuationDOF(2) = 0.00011111;%0.10628;
percentActuationDOF(3) = 7.8565e-05 ;%0.076607;
steps = 5;
[kinematicsData] = SequentialKinematics(flatFoldableTruss, steps, ...
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
anglesDeg
R = rotz(rad2deg(atan2(coords(end, 2),coords(end, 1))));
coords = [coords, zeros(size(coords,1), 1)] * R;
coords = coords(:, 1:2);
printMatrix2Col(coords);
function printMatrix2Col(A)
%PRINTMATRIX2COL Print an n-by-2 matrix in MATLAB row format.
%
% Example output:
% [0, 0;
%  0.5, 0.866025403784439;
%  1, 0]

    if ~ismatrix(A) || size(A,2) ~= 2
        error('Input must be an n-by-2 matrix.');
    end

    fprintf('[');
    for i = 1:size(A,1)
        if i == 1
            fprintf('%g, %g', A(i,1), A(i,2));
        else
            fprintf(';\n %g, %g', A(i,1), A(i,2));
        end
    end
    fprintf('];\n');
end