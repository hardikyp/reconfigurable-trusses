clc; clear; close all;
addpath("SequentialKinematics\");
addpath("TopologyOptimization\");
addpath("TraditionalTrusses\");
addpath("ModifiedTopOptTruss\");
addpath("Videos\");

trussName = "L-shape Beam";
nodes = [0, 1;
         0, 2;
         0.5, 1.5;
         1, 0;
         1, 1/sqrt(2);
         1, 1;
         1, 2;
         2, 0;
         2, 1;
         3, 1];

connectivity = [1, 2;
                1, 3;
                1, 4;
                1, 6;
                2, 7;
                3, 7;
                4, 5;
                4, 8;
                5, 6;
                6, 7;
                6, 8;
                6, 9;
                8, 10;
                9, 10];

supports = [2, 1, 1;
            7, 1, 1];

externalForces = [10, 0, -1];

barAreas = ones(size(connectivity, 1), 1);

% Member force is a dummy veriable for identifying tensile members
memberForces = [-1; 1; -1; 1; 1; 1; 1; -1; 1; 1; 1; 1; -1; 1];

D = [nodes(connectivity(:, 2), 1) - nodes(connectivity(:, 1), 1), ...
     nodes(connectivity(:, 2), 2) - nodes(connectivity(:, 1), 2)];
barLengths = sqrt(D(:, 1).^2 + D(:, 2).^2);
D = [D(:, 1)./barLengths, D(:, 2)./barLengths];

newNodeHist = [3, 1, 7, 1, 7, 2, 7;
               5, 4, 6, 3, 6, 1, 6;
               9, 6, 10, 5, 6, 6, 8];

triangles = [1, 2, 7;
             1, 6, 7;
             1, 4, 6;
             4, 6, 8;
             6, 8, 10];

flatFoldableTruss.nodes = nodes;
flatFoldableTruss.connectivity = connectivity;
flatFoldableTruss.supports = supports;
flatFoldableTruss.externalForces = externalForces;
flatFoldableTruss.barAreas = barAreas;
flatFoldableTruss.memberForces = memberForces;
flatFoldableTruss.D = D;
flatFoldableTruss.barLengths = barLengths;
flatFoldableTruss.newNodeNum = newNodeHist(:, 1);
flatFoldableTruss.newNodeHist = newNodeHist;
flatFoldableTruss.triangles = triangles;

numDOF = length(flatFoldableTruss.newNodeNum);
percentActuationDOF = ones(1, numDOF) * 0.85;
percentActuationDOF(end) = -percentActuationDOF(end);
steps = 30;
[kinematicsData] = SequentialKinematics(flatFoldableTruss, steps, ...
                                        percentActuationDOF);
% PlotSequentialKinematics(kinematicsData, trussName);
PlotConvexHull(kinematicsData, trussName, numDOF, steps);