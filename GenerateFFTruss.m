function [flatFoldableTruss] = GenerateFFTruss(processedData)
numTriangles = processedData.numTriangles;
triangles = processedData.triangles;
nodes = processedData.nodes;
memberForces = processedData.memberForces;
connectivity = processedData.connectivity;
barAreas = processedData.barAreas;
supports = processedData.supports;
externalForces = processedData.externalForces;
D = processedData.D;
newNodeHist = [0, 0, 0, 0, 0, 0, 0];

for triangle = 1:numTriangles
    [nodeI, nodeJ, nodeK] = deal(triangles(triangle, 1), ...
                                 triangles(triangle, 2), ...
                                 triangles(triangle, 3));
    currentTriangle = [nodeI, nodeJ, nodeK];
    bars = [nodeI nodeJ;
            nodeI nodeK;
            nodeJ nodeK];

    if ~any(ismember(bars, newNodeHist(:, 2:3), "rows"))
        d = [nodes(bars(:, 2), 1) - nodes(bars(:, 1), 1), ...
             nodes(bars(:, 2), 2) - nodes(bars(:, 1), 2)];

        triangleSides = sqrt(d(:, 1).^2 + d(:, 2).^2);
        triangleMemberForces = memberForces(ismember(connectivity, bars, 'rows'));
        triangleMemberAreas = barAreas(ismember(connectivity, bars, 'rows'));

        [triangleSides, sortIdx] = sort(triangleSides);
        triangleMemberForces = triangleMemberForces(sortIdx);
        triangleMemberAreas = triangleMemberAreas(sortIdx);
        bars = bars(sortIdx, :);

        % Initialize bar and triangle attributes
        previousTriangle = 0;
        nextTriangle = 0;
        groundBar = 0;
        avoidBar = 0;
        groundIdx = 0;
        avoidIdx = 0;

        % Find if structural support is present in triangle
        if nnz(ismember(currentTriangle, supports(:, 1))) > 1
            groundBar = sort(intersect(currentTriangle, supports(:, 1)))';
            groundIdx = find(ismember(bars, groundBar, 'rows'));
        end

        if numTriangles > 1
            if triangle == 1
                nextTriangle = triangles(triangle + 1, :);
                % If two nodes of the current triangle are not 
                % support nodes
                if (groundBar == 0) & (groundIdx == 0)
                    % Ground bar needs to be the bar that shares 
                    % a node with one of the supports & is compressive
                    suppNode = intersect(currentTriangle, supports(:, 1));
                    [rowIdx, ~] = find(bars==suppNode);
                    groundIdx = intersect(rowIdx, find(triangleMemberForces < 0));
                    groundBar = bars(groundIdx, :);
                end
                % Check if the next triangle shares a side with the
                % current triangle
                if ~isscalar(intersect(currentTriangle, nextTriangle))
                    avoidBar = sort(intersect(currentTriangle, nextTriangle));
                    avoidIdx = find(ismember(bars, avoidBar, 'rows'));
                end
            elseif triangle == numTriangles
                previousTriangle = triangles(triangle - 1, :);
                if (groundBar == 0) & (groundIdx == 0)
                    % Check if the previous triangle shares a side with the
                    % current triangle
                    if ~isscalar(intersect(currentTriangle, previousTriangle))
                        groundBar = sort(intersect(currentTriangle, previousTriangle));
                        groundIdx = find(ismember(bars, groundBar, 'rows'));
                    end
                end
            else
                previousTriangle = triangles(triangle - 1, :);
                if (groundBar == 0) & (groundIdx == 0)
                    %%%%% ADD A CONDITION TO MAKE GROUND BAR
                    %%%%% ANY BAR CONNECTED TO THE STRUCTURE
                    %%%%% EG. CANTILEVER 3,2,2,1
                    if ~isscalar(intersect(currentTriangle, previousTriangle))
                        groundBar = sort(intersect(currentTriangle, previousTriangle));
                        groundIdx = find(ismember(bars, groundBar, 'rows'));
                    end
                end
                nextTriangle = triangles(triangle + 1, :);
                % Check if the next triangle shares a side with the
                % current triangle
                if ~isscalar(intersect(currentTriangle, nextTriangle))
                    avoidBar = sort(intersect(currentTriangle, nextTriangle));
                    avoidIdx = find(ismember(bars, avoidBar, 'rows'));
                end
            end
        end

        candidateIdx = find(triangleMemberForces > 0);
        candidateIdx = setdiff(candidateIdx, [groundIdx, avoidIdx]);
        
        if ~isscalar(candidateIdx)
            [~, ind] = max(triangleMemberAreas(candidateIdx));
            candidateIdx = candidateIdx(ind);
        end

        if triangle < numTriangles
            [nextTriangleSimilar, nextRatio] = CheckSimilarTriangles( ...
                    currentTriangle, ...
                    nextTriangle, ...
                    nodes);
            if (nextTriangleSimilar) && (nextRatio - 1 < 1e-4) && (triangleMemberForces(avoidIdx) > 0)
                candidateIdx = avoidIdx;
            end
        end

        if length(unique(triangleSides)) == 3
            % Scalene triangle
            l1 = triangleSides(1);
            l2 = triangleSides(2);
            l3 = triangleSides(3);
            if candidateIdx == 1
                s = (l1 + l2 - l3) / 2;
            elseif candidateIdx == 2
                s = (l1 + l2 - l3) / 2;
            elseif candidateIdx == 3
                s = (l1 - l2 + l3) / 2;
            end

        elseif length(unique(triangleSides)) == 2
            % Isoceles triangle
            l1 = triangleSides(1);
            l2 = triangleSides(3);
            if triangleSides(2) == triangleSides(3)
                s = l1 / 2;
            elseif triangleSides(1) == triangleSides(2)
                if candidateIdx == 3
                    s = l2 / 2;
                else
                    s = l1 - l2 / 2;
                end
            end

        else
            % Equilateral triangle
            l1 = mean(triangleSides);
            s = l1 / 2;
        end
        
        if ~isempty(candidateIdx)
            idx = find(all(connectivity == bars(candidateIdx, :), 2));
            
            if any(groundBar)
                commonNode = intersect(bars(candidateIdx, :), groundBar);
                if commonNode == connectivity(idx, 2)
                    newNodes = nodes(commonNode, :) - s * D(idx, :);
                else
                    newNodes = nodes(commonNode, :) + s * D(idx, :);
                end
            else
                newNodes = nodes(bars(candidateIdx, 1), :) + s * D(idx, :);
            end
            
            idxNewNode = find((nodes(:, 1) < newNodes(1)) == 0, 1, "first");
            if (newNodes(1) == nodes(idxNewNode, 1))
                if newNodes(2) < nodes(idxNewNode, 2)
                    nodeNum = idxNewNode;
                else
                    nodeNum = idxNewNode + 1;
                end
            else
                nodeNum = idxNewNode;
            end
            % nodeNum = bars(candidateIdx, 2);

            % Add new node in nodes
            nodes = cat(1, ...
                        nodes(1:(nodeNum - 1), :), ...
                        newNodes, ...
                        nodes(nodeNum:end, :));
            % % Add new node in nodes
            % nodes = cat(1, ...
            %             nodes(1:(bars(candidateIdx, 2) - 1), :), ...
            %             newNodes, ...
            %             nodes(bars(candidateIdx, 2):end, :));

            if nodeNum > bars(groundIdx, 2)
                newNodeHist = [newNodeHist;
                               nodeNum, bars(candidateIdx, 1), ...
                               bars(candidateIdx, 2) + 1, triangle, ...
                               commonNode, bars(groundIdx, :)];
            else
                newNodeHist = [newNodeHist;
                               nodeNum, bars(candidateIdx, 1), ...
                               bars(candidateIdx, 2) + 1, triangle, ...
                               commonNode, bars(groundIdx, 1), ...
                               bars(groundIdx, 2) + 1];
            end
            rmvIdx = find(ismember(connectivity, bars(candidateIdx, :), "rows"));
            connectivity(rmvIdx, :) = [];
            
            candidateArea = barAreas(idx);
            barAreas(rmvIdx) = [];

            candidateForce = memberForces(idx);
            memberForces(rmvIdx) = [];

            [I, J] = find(connectivity >= nodeNum);
            for i = 1:length(I)
                if ~isequal(connectivity(I(i), :), bars(candidateIdx, :))
                    connectivity(I(i), J(i)) = connectivity(I(i), J(i)) + 1;
                end
            end


            % Modify connectivity of old system
            % [I, J] = find(connectivity > nodeNum);
            % for i = 1:length(I)
            %     if ~isequal(connectivity(I(i), :), bars(candidateIdx, :))
            %         connectivity(I(i), J(i)) = connectivity(I(i), J(i)) + 1;
            %     end
            % end


            % connectivity = cat(1, ...
            %                    connectivity, ...
            %                    [nodeNum, bars(candidateIdx, 2) + 1]);
            connectivity = cat(1, ...
                               connectivity, ...
                               [bars(candidateIdx, 1), nodeNum], ...
                               [nodeNum, bars(candidateIdx, 2) + 1]);
            [connectivity, conSort] = sortrows(connectivity);

            % Add bar area information
            % barAreas = [barAreas; 
            %             barAreas(idx) * ones(size(nodeNum, 1), 1)];
            barAreas = [barAreas; 
                        candidateArea * ones(2, 1)];
            barAreas = barAreas(conSort);

            % Add member force information
            % memberForces = [memberForces;
            %                 memberForces(idx) * ones(size(nodeNum, 1), 1)];
            memberForces = [memberForces;
                            candidateForce * ones(2, 1)];
            memberForces = memberForces(conSort);

            % Change triangle information
            [P, Q] = find(triangles >= nodeNum(1));
            for i = 1:length(P)
                triangles(P(i), Q(i)) = triangles(P(i), Q(i)) + 1;
            end

            % Find bar lengths and directional cosigns
            D = [nodes(connectivity(:, 2), 1) - nodes(connectivity(:, 1), 1), ...
                nodes(connectivity(:, 2), 2) - nodes(connectivity(:, 1), 2)];
            barLengths = sqrt(D(:, 1).^2 + D(:, 2).^2);
            D = [D(:, 1)./barLengths, D(:, 2)./barLengths];

            % Store external force information
            idxF = find(externalForces(:, 1) >= nodeNum(1));
            externalForces(idxF, 1) = externalForces(idxF, 1) + 1;

            % Remove unused support support nodes and store support information
            idxS = find(supports(:, 1) >= nodeNum(1));
            supports(idxS, 1) = supports(idxS, 1) + 1;
        end
    end
end

flatFoldableTruss.nodes = nodes;
flatFoldableTruss.barAreas = barAreas;
flatFoldableTruss.memberForces = memberForces;
flatFoldableTruss.connectivity = connectivity;
flatFoldableTruss.barLengths = barLengths;
flatFoldableTruss.supports = supports;
flatFoldableTruss.externalForces = externalForces;
flatFoldableTruss.D = D;
newNodeHist(~any(newNodeHist, 2), :) = [];
flatFoldableTruss.newNodeNum = newNodeHist(:, 1);
flatFoldableTruss.newNodeHist = newNodeHist;
flatFoldableTruss.triangles = triangles;
end