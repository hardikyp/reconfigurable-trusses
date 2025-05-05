function [] = PlotModifiedTruss(flatFoldableTruss)
nodes = flatFoldableTruss.nodes;
connectivity = flatFoldableTruss.connectivity;
barAreas = flatFoldableTruss.barAreas;
memberForces = flatFoldableTruss.memberForces;
externalForces = flatFoldableTruss.externalForces; 
supports = flatFoldableTruss.supports;
newNodeNum = flatFoldableTruss.newNodeNum;

fprintf('-PLOT- Modified flat foldable truss structure plotted\n');
figure('Name', 'Flat Foldable Truss', 'NumberTitle', 'off')
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
for n = 1:size(nodes, 1)
    if ismember(n, newNodeNum)
        plot(nodes(n, 1), nodes(n, 2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    else
        plot(nodes(n, 1), nodes(n, 2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    end
end

% -- Plot arrow for applied load -- %
arrow([nodes(externalForces(:, 1), 1), nodes(externalForces(:, 1), 2)], ...
      [(nodes(externalForces(:, 1), 1) + externalForces(:, 2) * max(nodes(:, 1)) / 10), ... 
       (nodes(externalForces(:, 1), 2) + externalForces(:, 3) * max(nodes(:, 2)) / 10)], ...
      'Color', 'r');

% -- Plot support conditions -- %
for i = 1:size(supports, 1)
    if all(supports(i, 2:3))
        plot(nodes(supports(i, 1), 1), nodes(supports(i, 1), 2), 'b>', ...
             'MarkerSize', 10, 'MarkerFaceColor', 'b');
    else
        plot(nodes(supports(i, 1), 1), nodes(supports(i, 1), 2), 'bo', ...
             'MarkerSize', 10, 'MarkerFaceColor', 'b');
    end
end
zoom out
end