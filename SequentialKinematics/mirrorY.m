function mirroredCoords = mirrorY(coords, anchorPt)
% Mirror about Y axis
coords(:, 1) = -1 * coords(:, 1);
% Translate relative to the anchor point
if length(anchorPt) == 3
    mirroredCoords = coords + 2 * [anchorPt(1), 0, 0];
else
    mirroredCoords = coords + 2 * [anchorPt(1), 0];
end
end

