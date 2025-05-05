function [] = PlotConvexHull(kinematicsData, trussName, numDOF, steps)
fname = sprintf("Convex Hull %s", trussName);
displayText = sprintf("-PLOT- 2D Convex Hull of reconfigurable %s\n", trussName);
fprintf(displayText);
figure('Name', fname, 'NumberTitle', 'off');
nodeLoc = kinematicsData.kinematicsNodeLoc;
chArea = zeros(size(nodeLoc, 3), 1);
longestLength = zeros(size(nodeLoc, 3), 1);

for incr = 1:size(nodeLoc, 3)
    [chIdx, chArea(incr)] = convhull(nodeLoc(:, :, incr));
    d = pdist(nodeLoc(chIdx, :, incr));
    longestLength(incr) = max(d);
end

yyaxis left
plot(1:size(nodeLoc, 3), (chArea ./ chArea(1)), '-', LineWidth = 2);
ylabel("Area");
ylim([0 1.4])


yyaxis right
plot(1:size(nodeLoc, 3), (longestLength ./ longestLength(1)), '-', LineWidth = 2);
ylim([0 1.4])
ylabel("Longest length");

xlabel("Steps");
% xticks([0 5 10]);
% xticklabels({'x = 0', 'x = 5', 'x = 10'});

for i = 1:numDOF
    xline(i * (steps + 1));
end
end