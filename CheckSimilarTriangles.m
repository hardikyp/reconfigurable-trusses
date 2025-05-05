function [isSimilar, ratio] = CheckSimilarTriangles(triangle1, triangle2, nodes)
coordTriangle1 = nodes(triangle1, :);
coordTriangle2 = nodes(triangle2, :);

% Compute distances between each pair of vertices for both triangles
distances1 = pdist(coordTriangle1);
distances2 = pdist(coordTriangle2);

% Sort the distances in ascending order for each triangle
sorted_distances1 = sort(distances1);
sorted_distances2 = sort(distances2);

[~, ratio1] = convhull(coordTriangle1);
[~, ratio2] = convhull(coordTriangle2);

if cross(sorted_distances1, sorted_distances2) < 1e-5
    % The triangles are similar
    isSimilar = true;
    ratio = ratio2 / ratio1;
else
    % The triangles are not similar
    isSimilar = false;
    ratio = 0;
end

end
