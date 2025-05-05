function []= PlotGroundStructure(NODE, BARS, A, Cutoff, Ng, LOAD, SUPP)

figure('Name', 'GRAND v1.0 -- Zegard T, Paulino GH', 'NumberTitle', 'off')
hold on; axis equal; axis off; color = jet(Ng); axis('manual');
A = A / max(A); % Normalize to [0,1] areas
ind = find(A > Cutoff);
MyGroup = ceil(Ng * A(ind)); % Round up to the closest group of bars
Groups = cell(Ng, 1);       % Store the indices of similar bars
for i = 1:Ng
    Groups{i} = ind(MyGroup == i); 
end

for i = Ng:-1:1 % Plot each group of similar bars in a single plot call
    if ~isempty(Groups{i})
        XY = [NODE(BARS(Groups{i}, 1), :) NODE(BARS(Groups{i}, 2), :)];
        GroupArea = mean(A(Groups{i})); % Mean area for this group
        plot(XY(:, [1 3])', XY(:, [2 4])', 'LineWidth', ...
             5 * sqrt(GroupArea), 'Color', color(i, :))
    end
end

arrow([NODE(LOAD(:, 1), 1), NODE(LOAD(:, 1), 2)], ...
      [(NODE(LOAD(:, 1), 1) + LOAD(:, 2) * max(NODE(:, 1)) / 5), ... 
       (NODE(LOAD(:, 1), 2) + LOAD(:, 3) * max(NODE(:, 2)) / 5)], ...
      'Color', 'r');

plot(NODE(SUPP(:, 1), 1), NODE(SUPP(:, 1), 2), 'b>', 'MarkerSize', 8, ...
     'MarkerFaceColor', 'b');

% -- Plot nodes -- %
plot(NODE(:, 1), NODE(:, 2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

fprintf('-PLOT- Cutoff %g, Groups %g, Bars plotted %g\n', Cutoff, Ng, ...
        length(ind))