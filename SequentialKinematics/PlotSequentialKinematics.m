function [] = PlotSequentialKinematics(kinematicsData, trussName)
fprintf('-PLOT- Sequential kinematics of reconfigurable truss\n');
figure('Name', 'Sequential kinematics', 'NumberTitle', 'off');
videoFile = sprintf("Videos\\SequentialKinematics_%s.mp4", trussName);
v = VideoWriter(videoFile);
v.FrameRate = 15;
open(v);
nodeLoc = kinematicsData.kinematicsNodeLoc;
connectivity = kinematicsData.connectivity;
newNodeNum = kinematicsData.newNodeNum;
memberForces = kinematicsData.memberForces;
supports = kinematicsData.supports;
barAreas = kinematicsData.barAreas;
steps = [1:1:size(nodeLoc, 3), size(nodeLoc, 3):-1:1];
xMin = min(min(nodeLoc(:, 1, :))) - 1;
yMin = min(min(nodeLoc(:, 2, :))) - 1;
xMax = max(max(nodeLoc(:, 1, :))) + 1;
yMax = max(max(nodeLoc(:, 2, :))) + 1;
for incr = steps
    hold on; axis equal; axis off; color = [0 0.4470 0.7410; 0 0 0; 0.8500 0.3250 0.0980];
    xlim([xMin xMax]); ylim([yMin yMax]);
    % -- Plot bars -- %
    for bar = 1:size(connectivity, 1)
        if memberForces(bar) > 0
            plot([nodeLoc(connectivity(bar, 1), 1, incr), nodeLoc(connectivity(bar, 2), 1, incr)], ...
                 [nodeLoc(connectivity(bar, 1), 2, incr), nodeLoc(connectivity(bar, 2), 2, incr)], ...
                 'Color', color(1, :), 'LineWidth', 3 * barAreas(bar));
        elseif memberForces(bar) == 0
            plot([nodeLoc(connectivity(bar, 1), 1, incr), nodeLoc(connectivity(bar, 2), 1, incr)], ...
                 [nodeLoc(connectivity(bar, 1), 2, incr), nodeLoc(connectivity(bar, 2), 2, incr)], ...
                 'Color', color(2, :), 'LineWidth', 3 * barAreas(bar));
        else
            plot([nodeLoc(connectivity(bar, 1), 1, incr), nodeLoc(connectivity(bar, 2), 1, incr)], ...
                 [nodeLoc(connectivity(bar, 1), 2, incr), nodeLoc(connectivity(bar, 2), 2, incr)], ...
                 'Color', color(3, :), 'LineWidth', 3 * barAreas(bar));
        end
    end

    % -- Plot nodes -- %
    for n = 1:size(nodeLoc, 1)
        if ismember(n, newNodeNum)
            plot(nodeLoc(n, 1, incr), nodeLoc(n, 2, incr), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        else
            plot(nodeLoc(n, 1, incr), nodeLoc(n, 2, incr), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
        end
    end

    % -- Plot support conditions -- %
    for i = 1:size(supports, 1)
        if all(supports(i, 2:3))
            plot(nodeLoc(supports(i, 1), 1, incr), nodeLoc(supports(i, 1), 2, incr), 'b>', ...
                 'MarkerSize', 10, 'MarkerFaceColor', 'b');
        else
            plot(nodeLoc(supports(i, 1), 1, incr), nodeLoc(supports(i, 1), 2, incr), 'bo', ...
                 'MarkerSize', 10, 'MarkerFaceColor', 'b');
        end
    end

    % -- Write frame and clear for next frame -- %
    frame = getframe(gcf);
    writeVideo(v, frame);
    clf;
end
end