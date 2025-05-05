function [processedData] = FindTriangles(processedData)
    triangles = [0, 0, 0];
    triangleCount = 0;
    nodes = processedData.nodes;
    connectivity = processedData.connectivity;
    for node = 1:size(nodes, 1)
        [I, ~] = find(connectivity == node);
        for iBar = 1:length(I)
            nodeI = node;
            nodeJ = setdiff(connectivity(I(iBar), :), nodeI);
    
            % NodeI is connected to which nodes?
            nodeIConnectivity = setdiff(unique(connectivity(I, :)), ...
                                        [nodeI nodeJ]);
            
            % NodeJ is connected to which nodes?
            [P, ~] = find(connectivity == nodeJ);
            nodeJConnectivity = setdiff(unique(connectivity(P, :)), ...
                                        [nodeI nodeJ]);
            
            commonNodes = intersect(nodeIConnectivity, nodeJConnectivity);
    
            for i = 1:length(commonNodes)
                triangle = sort([nodeI, nodeJ, commonNodes(i)]);
                if ~ismember(triangle, triangles, 'rows')
                    triangleCount = triangleCount + 1;
                    triangles(triangleCount, :) = triangle;
                end
            end
        end
    end
    processedData.numTriangles = triangleCount;
    processedData.triangles = triangles;
end