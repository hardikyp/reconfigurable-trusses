%% Load Warren Truss Data
nodes = [0, 0;
         1, 1;
         2, 0;
         3, 1;
         4, 0;
         5, 1;
         6, 0];

connectivity = [1 ,2;
                1, 3;
                2, 3;
                2, 4;
                3, 4;
                3, 5;
                4, 5;
                4, 6;
                5, 6;
                5, 7;
                6, 7];

supports = [1, 1, 1;
            7, 1, 1];

externalForces = [1, 0, -1;
                  3, 0, -1;
                  5, 0, -1;
                  7, 0, -1];

barAreas = ones(size(connectivity, 1), 1);

% Member force is a dummy veriable for identifying tensile members
memberForces = [-1; 1; 1; -1; -1; 1; -1; -1; 1; 1; -1];

D = [nodes(connectivity(:, 2), 1) - nodes(connectivity(:, 1), 1), ...
     nodes(connectivity(:, 2), 2) - nodes(connectivity(:, 1), 2)];
barLengths = sqrt(D(:, 1).^2 + D(:, 2).^2);
D = [D(:, 1)./barLengths, D(:, 2)./barLengths];

inputData.nodes = nodes;
inputData.connectivity = connectivity;
inputData.supports = supports;
inputData.externalForces = externalForces;
inputData.barAreas = barAreas;
inputData.memberForces = memberForces;
inputData.D = D;