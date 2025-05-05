function [processedData] = ReducedOptimalGS(NODE, ...
                                            BARS, ...
                                            LOAD, ...
                                            SUPP, ...
                                            A, ...
                                            N, ...
                                            Nn, ...
                                            ColTol, ...
                                            Cutoff)
    %% Remove zero area bars, remove unused nodes and renumber nodes
    idx = A > Cutoff;
    barAreas = A(idx);
    connectivity = BARS(idx, :);
    memberForces = N(idx, :);
    uniqueNodes = unique(connectivity);
    nodes = NODE(uniqueNodes, :);
    nodesPurged = setdiff(1:Nn, uniqueNodes);
    
    while ~isempty(nodesPurged)
        [I, J] = find(connectivity >= nodesPurged(1));
        for i = 1:length(I)
            connectivity(I(i), J(i)) = connectivity(I(i), J(i)) - 1;
        end
        nodesPurged(1) = [];
        nodesPurged = nodesPurged - 1;
    end
    
    D = nodes(connectivity(:, 2), :) - nodes(connectivity(:, 1), :);
    barLengths = sqrt(sum(D.^2, 2));
    D = D ./ barLengths;
    
    %% Find nodes with collinear bars and swap out with single bar
    nodesToPurge = zeros(size(nodes, 1), 1);
    for n = 1:size(nodes, 1)
        [I, ~] = find(connectivity == n);
        if size(I, 1) == 2 && dot(D(I(1), :), D(I(2), :)) > ColTol
            nodesToPurge(n) = n;
        end
    end
    nodesToPurge = nonzeros(nodesToPurge);
    
    while ~isempty(nodesToPurge)
        nodeF = nodesToPurge(1);
        nodes(nodeF, :) = [];
        [I, J] = find(connectivity == nodeF);
        barAreas(I) = max(barAreas(I));
    
        barAreas(I(1)) = [];
        memberForces(I(1)) = [];
        connectivity(I(2), J(2)) = connectivity(I(1), 2);
        connectivity(I(1), :) = [];
        % reduce node number for subsequent nodes
        [P, Q] = find(connectivity >= nodeF);
        for i = 1:length(P)
            connectivity(P(i), Q(i)) = connectivity(P(i), Q(i)) - 1;
        end
        nodesToPurge(1) = [];
        nodesToPurge = nodesToPurge - 1;
    end

    %% Store external force information
    idxNode = find(ismember(nodes, NODE(LOAD(:, 1), :), 'rows'));
    idxForce = find(ismember(NODE(LOAD(:, 1), :), nodes, 'rows'));
    externalForces = [idxNode LOAD(idxForce, 2) LOAD(idxForce, 3)];
    
    %% Remove unused support support nodes and store support information
    idxNode = find(ismember(nodes, NODE(SUPP(:, 1), :), 'rows'));
    idxSupp = find(ismember(NODE(SUPP(:, 1), :), nodes, 'rows'));
    supports = [idxNode SUPP(idxSupp, 2) SUPP(idxSupp, 3)];

    %% Add zero force members between supports
    tempConnectivity = connectivity;
    tempBarAreas = barAreas;
    tempMemberForces = memberForces;
    suppConn = [];
    for i = 1:size(supports, 1) - 1
        N1 = supports(i, 1); N2 = supports(i + 1, 1);
        % Check if N1 and N2 are directly connected to each other
        if ~any(ismember([N1, N2], tempConnectivity, 'rows'))
            % Check for additional nodes on the line connecting
            % N1 and N2
            dN1 = nodes - nodes(N1, :);
            lN1 = sqrt(dN1(:, 1).^2 + dN1(:, 2).^2);
            dN1 = dN1 ./ lN1;
            collinearNodes = setdiff(find(ismember(dN1, dN1(N2, :), 'rows')), N2);
            % If nodes present between N1 and N2, connect the ones 
            % not already connected
            if ~isempty(collinearNodes)
                tempBar = reshape([N1, repelem(collinearNodes', 2), N2], 2, [])';
                tempBar = setdiff(tempBar, connectivity, 'rows');
            else
                tempBar = [N1, N2];
            end
            suppConn = [suppConn; tempBar];
        end
        % suppConn = [suppConn; tempBar];
    end
    suppConn = setdiff(suppConn, tempConnectivity, 'rows');
    tempConnectivity = [tempConnectivity; suppConn];
    tempBarAreas = [tempBarAreas; 
                    min(barAreas) * ones(size(suppConn, 1), 1)];
    tempMemberForces = [tempMemberForces;
                        zeros(size(suppConn, 1), 1)];
    [connectivity, sortIdx] = sortrows(tempConnectivity);
    barAreas = tempBarAreas(sortIdx);
    memberForces = tempMemberForces(sortIdx);

    %% Check and connect stray bar at ext. force node to nearest node
    tempConnectivity = connectivity;
    tempBarAreas = barAreas;
    tempMemberForces = memberForces;

    for i = 1:size(externalForces, 1)
        [I, ~] = find(connectivity == externalForces(i, 1));
        if isscalar(I)
            nodeF = externalForces(i, 1);
            existingConnectivity = setdiff(connectivity(I, :), nodeF);

            tempD = nodes - nodes(nodeF, :);
            tempDist = sqrt(tempD(:, 1).^2 + tempD(:, 2).^2);
            [sortedTempDist, sortIdx] = sort(tempDist, 1);
            uniqueDist = unique(sortedTempDist);
            closestNodes = sortIdx(sortedTempDist == uniqueDist(2));

            % Need to implement a check on crossing bars for a
            % robust way of reconnecting stray bar in future
            if isscalar(closestNodes)
                closestNodes = sortIdx(sortedTempDist == uniqueDist(3));
            end

            for n = 1:length(closestNodes)
               if closestNodes(n) ~= existingConnectivity
                   tempConnectivity = [tempConnectivity; ...
                                       sort([nodeF, closestNodes(n)])];
                   tempBarAreas = [tempBarAreas; min(barAreas)];
                   tempMemberForces = [tempMemberForces; 0];
               end
            end
        end
    end
    [connectivity, sortIdx] = sortrows(tempConnectivity);
    barAreas = tempBarAreas(sortIdx);
    memberForces = tempMemberForces(sortIdx);

    %% Recalc. bar lengths and directional unit vectors
    D = nodes(connectivity(:, 2), :) - nodes(connectivity(:, 1), :);
    barLengths = sqrt(sum(D.^2, 2));
    D = D ./ barLengths;

    %% Store data in a struct
    processedData.barAreas = barAreas;
    processedData.connectivity = connectivity;
    processedData.memberForces = memberForces;
    processedData.nodes = nodes;
    processedData.barLengths = barLengths;
    processedData.D = D;
    processedData.externalForces = externalForces;
    processedData.supports = supports;
end