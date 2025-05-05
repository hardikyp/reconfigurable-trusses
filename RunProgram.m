clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath("SequentialKinematics\");
addpath("TopologyOptimization\");
addpath("TraditionalTrusses\");
addpath("ModifiedTopOptTruss\");
addpath("Videos\");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trussType = "traditional";
% trussType = "topologyOptimized";

trussName = "WarrenTruss";
% trussName = "FinkTruss";
% trussName = "KingPostTruss";
% trussName = "ScissorsTruss";

% Modified topology optimized trusses
% trussName = "SerpentineTrussN";

if trussType == "traditional"
    filename = trussName + ".m";
    run(filename);

elseif trussType == "topologyOptimized"
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRAND - Ground Structure Analysis and Design Code.
    % Tomas Zegard, Glaucio H Paulino - Version 1.0, Dec-2013
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% === MESH GENERATION LOADS/BCS ======================================
    kappa = 1.0; ColTol = 0.999999;
    Cutoff = 0.002; Ng = 50; % Plot: Member Cutoff & Number of plot groups
    
    % ------------- STRUCTURED-ORTHOGONAL MESH GENERATION -----------------
    % trussName = 'Bridge';
    % trussName = 'Cantilever';
    % [NODE, ELEM, SUPP, LOAD] = StructDomain(4, 1, 4, 1, trussName);
    % Lvl = 1; RestrictDomain = []; % No restriction for box domain

    % ------------- ORTHOGONAL L MESH GENERATION --------------------------
    % trussName = 'L-truss';
    % [NODE, ELEM, SUPP, LOAD] = LMesh(2, 2, 2, 2);
    % Lvl = 2; RestrictDomain = @RestrictL;

    % ------------- ORTHOGONAL SERPENTINE MESH GENERATION -----------------
    % [NODE, ELEM, SUPP, LOAD] = SerpentineMesh(4, 1, 8, 1);
    % Lvl = 1; RestrictDomain = @RestrictSerpentine;
    
    %% === GROUND STRUCTURE METHOD ========================================
    % PlotPolyMesh(NODE, ELEM, SUPP, LOAD);
    [BARS] = GenerateGS(NODE, ELEM, Lvl, RestrictDomain, ColTol);   
    Nn = size(NODE, 1); Ne = length(ELEM); Nb = size(BARS, 1);
    [BC] = GetSupports(SUPP);                                       
    [BT, L] = GetMatrixBT(NODE, BARS, BC, Nn, Nb);                  
    [F] = GetVectorF(LOAD, BC, Nn);
    
    fprintf('Mesh: Elements %d, Nodes %d, Bars %d, Level %d\n', ...
            Ne, Nn, Nb, Lvl)
    BTBT = [BT -BT]; LL = [L; kappa * L]; sizeBTBT = whos('BTBT');
    clear BT L;
    fprintf('Matrix [BT -BT]: %d x %d in %gMB (%gGB full)\n', ...
            length(F), length(LL), sizeBTBT.bytes/2^20, ...
            16 * (2 * Nn) * Nb / 2^30)
    
    tic, [S, vol, exitflag] = linprog(LL, [], [], BTBT, F, zeros(2 * Nb, 1));
    fprintf('Objective V = %f\nlinprog CPU time = %g s\n', vol, toc);
    
    S = reshape(S, [], 2);            % Separate slack variables
    A = S(:, 1) + kappa * S(:, 2);    % Get cross-sectional areas
    N = S(:, 1) - S(:, 2);            % Get member forces
    
    %% === PLOTTING =======================================================
    PlotGroundStructure(NODE, BARS, A, Cutoff, Ng, LOAD, SUPP);
    PlotBoundary(ELEM, NODE);

    %% === OBTAIN REDUCED OPTIMAL STRUCTURE ===============================
    [inputData] = ReducedOptimalGS(NODE, BARS, LOAD, SUPP, A, N, Nn, ...
                                       ColTol, Cutoff);
    
    %% === PLOT REDUCED OPTIMAL STRUCTURE =================================
    PlotROS(ELEM, NODE, inputData);
end

%% === IDENTIFY TRIANGULAR ELEMENTS =======================================
[processedData] = FindTriangles(inputData);

%% === LOCATION FOR NEW NODES ENABLING FLAT FOLDABILITY ===================
[flatFoldableTruss] = GenerateFFTruss(processedData);

%% === PLOT MODIFIED TRUSS STRUCTURE ======================================
PlotModifiedTruss(flatFoldableTruss);

%% === KINEMATIC ANALYSIS =================================================
numDOF = length(flatFoldableTruss.newNodeNum);
percentActuationDOF = ones(1, numDOF) * 0.85;
if trussName == "ScissorsTruss"
    percentActuationDOF([1, 4]) = -percentActuationDOF([1, 4]);
elseif trussName == "FinkTruss"
    percentActuationDOF(end) = 0;
    elseif trussName == "SerpentineTruss"
    percentActuationDOF(end) = -percentActuationDOF(end);
end
steps = 30;
[kinematicsData] = SequentialKinematics(flatFoldableTruss, steps, ...
                                        percentActuationDOF);

%% === PLOT KINEMATICS ====================================================
% PlotSequentialKinematics(kinematicsData, trussName);
PlotConvexHull(kinematicsData, trussName, numDOF, steps);