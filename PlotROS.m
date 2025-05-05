function [] = PlotROS(ELEM, NODE, processedData)
nodes = processedData.nodes;
connectivity = processedData.connectivity;
barAreas = processedData.barAreas;
memberForces = processedData.memberForces;
externalForces = processedData.externalForces; 
supports = processedData.supports;

fprintf('-PLOT- Reduced optimal ground structure plotted\n');
figure('Name', 'Processed Ground Structure', 'NumberTitle', 'off')
hold on; axis equal; axis off; color = [0 0.4470 0.7410; 0 0 0; 0.8500 0.3250 0.0980]; axis('manual');

% -- Plot bars -- %
for bar = 1:size(connectivity, 1)
    if memberForces(bar) > 0
        plot([nodes(connectivity(bar, 1), 1), nodes(connectivity(bar, 2), 1)], ...
             [nodes(connectivity(bar, 1), 2), nodes(connectivity(bar, 2), 2)], ...
             'Color', color(1, :), 'LineWidth', 3 * barAreas(bar));
    elseif memberForces(bar) == 0
        plot([nodes(connectivity(bar, 1), 1), nodes(connectivity(bar, 2), 1)], ...
             [nodes(connectivity(bar, 1), 2), nodes(connectivity(bar, 2), 2)], ...
             'Color', color(2, :), 'LineWidth', 3 * barAreas(bar));
    else
        plot([nodes(connectivity(bar, 1), 1), nodes(connectivity(bar, 2), 1)], ...
             [nodes(connectivity(bar, 1), 2), nodes(connectivity(bar, 2), 2)], ...
             'Color', color(3, :), 'LineWidth', 3 * barAreas(bar));
    end
end
% -- Plot nodes -- %
plot(nodes(:, 1), nodes(:, 2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

% -- Plot arrow for applied load -- %
arrow([nodes(externalForces(:, 1), 1), nodes(externalForces(:, 1), 2)], ...
      [(nodes(externalForces(:, 1), 1) + externalForces(:, 2) * max(nodes(:, 1)) / 5), ... 
       (nodes(externalForces(:, 1), 2) + externalForces(:, 3) * max(nodes(:, 2)) / 5)], ...
      'Color', 'r');

% -- Plot support conditions -- %
for i = 1:size(supports, 1)
    if all(supports(i, 2:3))
        plot(nodes(supports(i, 1), 1), nodes(supports(i, 1), 2), 'b>', ...
             'MarkerSize', 8, 'MarkerFaceColor', 'b');
    else
        plot(nodes(supports(i, 1), 1), nodes(supports(i, 1), 2), 'bo', ...
             'MarkerSize', 8, 'MarkerFaceColor', 'b');
    end
end

% -- Plot boundary -- %
Nn = size(NODE, 1); Ne = length(ELEM); NpE = cellfun(@numel, ELEM);

FACE = sparse([], [], [], Nn, Nn, sum(NpE));
for i = 1:Ne
    MyFACE = [ELEM{i}; ELEM{i}(2:end) ELEM{i}(1)];
    for j = 1:NpE(i)
        if FACE(MyFACE(1, j), MyFACE(2,j)) == 0 % New edge - Flag it
            FACE(MyFACE(1, j), MyFACE(2, j)) = i;
            FACE(MyFACE(2, j), MyFACE(1, j)) = -i;
        elseif isnan(FACE(MyFACE(1, j), MyFACE(2, j)))
            error(sprintf('Edge [%d %d] found in >2 elements', MyFACE(:, j)))
        else % Edge belongs to 2 elements: inside domain. Lock it.
            FACE(MyFACE(1, j), MyFACE(2, j)) = NaN;
            FACE(MyFACE(2, j), MyFACE(1, j)) = NaN;
        end
    end
end
[BOUND(:, 1), BOUND(:, 2)] = find(FACE > 0);
BOUND(:, 3) = FACE(sub2ind(size(FACE), BOUND(:, 1), BOUND(:, 2)));
% plot([NODE(BOUND(:, 1), 1) NODE(BOUND(:, 2), 1)]',[NODE(BOUND(:, 1), 2) NODE(BOUND(:, 2), 2)]', '--k');
end