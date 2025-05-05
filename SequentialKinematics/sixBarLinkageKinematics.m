function [individualDOFNodeLoc, dPhi] = sixBarLinkageKinematics(linkageCoordinates, ...
                                                                 linkageConnectivity, ...
                                                                 groundLink, ...
                                                                 originNode, addedNode, ...
                                                                 pctActuation, steps)

lclOrigin = linkageCoordinates(originNode, :);
nonOrgGnd = setdiff(linkageConnectivity(groundLink, :), originNode);
d = [linkageCoordinates(nonOrgGnd, 1) - ...
     linkageCoordinates(originNode, 1), ...
     linkageCoordinates(nonOrgGnd, 2) - ...
     linkageCoordinates(originNode, 2)];
l = sqrt(sum(d.^2));
d = d ./ l;
angleGround = atan2(d(2), d(1));
% angleGround = mod(angleGround, 2*pi);

angleAddedNode = atan2(linkageCoordinates(addedNode, 2) - linkageCoordinates(originNode, 2), ...
                     linkageCoordinates(addedNode, 1) - linkageCoordinates(originNode, 1));
% angleAddedNode = mod(angleAddedNode, 2*pi);

if angleGround > 0
    if (angleAddedNode < angleGround) && (angleAddedNode > angleGround - pi) 
        lclAxes = rotz(rad2deg(pi + angleGround));
    else
        lclAxes = rotz(rad2deg(angleGround));
    end
else
    if (angleAddedNode <= angleGround + pi)
        angleGround = mod(angleGround, 2*pi);
        lclAxes = rotz(rad2deg(angleGround));
    else
        angleGround = mod(angleGround, 2*pi);
        lclAxes = rotz(rad2deg(angleGround + pi));
    end
end

% if angleGround <= pi
%     if (angleAddedNode > angleGround) && (angleAddedNode <= angleGround + pi) 
%         lclAxes = rotz(rad2deg(angleGround));
%     else
%         lclAxes = rotz(rad2deg(pi + angleGround));
%     end
% else
%     if (angleAddedNode > angleGround) && (angleAddedNode <= angleGround + pi)
%         lclAxes = rotz(rad2deg(angleGround));
%     else
%         lclAxes = rotz(rad2deg(angleGround - pi));
%     end
% end

lclLinkageCoordinates = global2localcoord(linkageCoordinates', 'rr', ...
                                          lclOrigin', lclAxes);
lclLinkageCoordinates = lclLinkageCoordinates(1:2, :)';
d = [lclLinkageCoordinates(nonOrgGnd, 1) - ...
     lclLinkageCoordinates(originNode, 1), ...
     lclLinkageCoordinates(nonOrgGnd, 2) - ...
     lclLinkageCoordinates(originNode, 2)];
l = sqrt(sum(d.^2));
d = d ./ l;
if abs(pi - abs(atan2(d(2), d(1)))) < 1e-6
    doFlip = true;
else
    doFlip = false;
end

if doFlip
    lclLinkageCoordinates = mirrorY(lclLinkageCoordinates, ...
                                    lclLinkageCoordinates(originNode, :));
end

% Find bar lengths and directional cosines for each bar
D = [lclLinkageCoordinates(linkageConnectivity(:, 2), 1) - ...
     lclLinkageCoordinates(linkageConnectivity(:, 1), 1), ...
     lclLinkageCoordinates(linkageConnectivity(:, 2), 2) - ...
     lclLinkageCoordinates(linkageConnectivity(:, 1), 2)];
linkageSides = sqrt(D(:, 1).^2 + D(:, 2).^2);
D = [D(:, 1)./linkageSides, D(:, 2)./linkageSides];

% a, b, c, d, x and y are the links of six bar linkage
% d is always the ground link
d = linkageSides(groundLink);

% c is the output link connected to groundLink
[I, ~] = find(linkageConnectivity == ...
              setdiff(linkageConnectivity(groundLink, :), originNode));
cIdx = setdiff(I, groundLink);
c = linkageSides(cIdx);

% a is the input link. Connects addedNode and originNode
aIdx = find(ismember(linkageConnectivity, ...
                     sort([originNode, addedNode]), "rows"));
a = linkageSides(aIdx);

% b is the coupler link connecting a & c.
nodeC = setdiff(linkageConnectivity(cIdx, :), nonOrgGnd);
bIdx = find(ismember(linkageConnectivity, sort([addedNode, nodeC]), "rows"));
b = linkageSides(bIdx);

% x is the input link of additional part
nodeE = setdiff(1:size(linkageCoordinates, 1), [originNode, nonOrgGnd, addedNode, nodeC]);
xIdx = find(ismember(linkageConnectivity, sort([originNode, nodeE]), "rows"));
x = linkageSides(xIdx);

% y is the coupler link of additional part
yIdx = find(ismember(linkageConnectivity, sort([nodeC, nodeE]), "rows"));
y = linkageSides(yIdx);



% Initialize linkage angles
% theta = [theta2, theta3, theta4, theta2Prime, theta3Prime]
theta = zeros(steps + 1, 5);
linkageNodeLoc = zeros(size(linkageCoordinates, 1), ...
    size(linkageCoordinates, 2), steps + 1);
linkageNodeLoc(:, :, 1) = linkageCoordinates;
theta2 = atan2(lclLinkageCoordinates(addedNode, 2) - lclLinkageCoordinates(originNode, 2), ...
               lclLinkageCoordinates(addedNode, 1) - lclLinkageCoordinates(originNode, 1));
theta3 = atan2(lclLinkageCoordinates(nodeC, 2) - lclLinkageCoordinates(addedNode, 2), ...
               lclLinkageCoordinates(nodeC, 1) - lclLinkageCoordinates(addedNode, 1));
theta4 = atan2(lclLinkageCoordinates(nodeC, 2) - lclLinkageCoordinates(nonOrgGnd, 2), ...
               lclLinkageCoordinates(nodeC, 1) - lclLinkageCoordinates(nonOrgGnd, 1));
theta5 = atan2(lclLinkageCoordinates(nodeE, 2) - lclLinkageCoordinates(originNode, 2), ...
               lclLinkageCoordinates(nodeE, 1) - lclLinkageCoordinates(originNode, 1));
theta6 = atan2(lclLinkageCoordinates(nodeE, 2) - lclLinkageCoordinates(nodeC, 2), ...
               lclLinkageCoordinates(nodeE, 1) - lclLinkageCoordinates(nodeC, 1));
% theta(1, :) = [atan2(D(aIdx, 2), D(aIdx, 1)), ...
%                atan2(D(bIdx, 2), D(bIdx, 1)), ...
%                atan2(D(cIdx, 2), D(cIdx, 1)), ...
%                atan2(D(xIdx, 2), D(xIdx, 1)), ...
%                atan2(D(yIdx, 2), D(yIdx, 1))];
theta(1, :) = [theta2, theta3, theta4, theta5, theta6];

if pctActuation > 0
    theta(:, 1) = linspace(theta(1, 1), (1 - pctActuation) * theta(1, 1), ...
                           steps + 1)';
elseif pctActuation < 0
    theta(:, 1) = linspace(theta(1, 1), -pctActuation * (pi), steps + 1)';
else
    % theta(:, 1) = linspace(theta(1, 1), -theta(1, 1), steps + 1)';
    theta(:, 1) = linspace(theta(1, 1), theta(1, 1) - 0.85 * (theta(1, 1) + pi), steps + 1)'; 
end

lclLinkageNodeLoc = zeros(5, 3);
initGuessE = [lclLinkageCoordinates(nodeE, 1), ...
              lclLinkageCoordinates(nodeE, 2)];

for inc = 2:steps+1
    % Find new theta3 and theta4
    r = d - a * cos(theta(inc, 1));
    s = a * sin(theta(inc, 1));
    cos_delta = (b^2 + c^2 - r^2 - s^2) / (2 * b * c);
    delta = acos(cos_delta);
    g = b - c * cos(delta);
    h = c * sin(delta);
    theta(inc, 2) = atan2((h * r - g * s), (g * r + h * s));
    theta(inc, 3) = theta(inc, 2) + delta;

    % Node A index (local)
    lclLinkageNodeLoc(originNode, :) = [0, 0, 0];
    % Node D index (local)
    idxD = setdiff(linkageConnectivity(groundLink, :), originNode);
    lclLinkageNodeLoc(idxD, :) = lclLinkageNodeLoc(originNode, :) + [d, 0, 0];
    % Node B index (local)
    idxB = setdiff(linkageConnectivity(aIdx, :), originNode);
    lclLinkageNodeLoc(idxB, :) = lclLinkageNodeLoc(originNode, :) + ...
        [a * cos(theta(inc, 1)), a * sin(theta(inc, 1)), 0];
    % Node C index (local)
    idxC = setdiff(linkageConnectivity(cIdx, :), idxD);
    lclLinkageNodeLoc(idxC, :) = lclLinkageNodeLoc(idxD, :) + ...
        [c * cos(theta(inc, 3)), c * sin(theta(inc, 3)), 0];
    % if ~(doFlip)
    %     rP = d - x * cos(pi - theta(inc, 3));
    %     sP = a * sin(pi - theta(inc, 3));
    %     deltaP = acos((y^2 + x^2 - rP^2 -sP^2) / (2 * x * y));
    %     gP = y - x * cos(deltaP);
    %     hP = c * sin(deltaP);
    %     t3P = atan2((hP * rP - gP * sP), (gP * rP + hP * sP));
    %     t4P = t3P + deltaP;
    %     theta(inc, 4) = pi - t4P;
    % end

    % Node E index (local)
    idxE = setdiff(linkageConnectivity(xIdx, :), originNode);
    fun = @(P) [(P(1) - lclLinkageNodeLoc(originNode, 1))^2 + ...
                (P(2) - lclLinkageNodeLoc(originNode, 2))^2 - x^2;  % dist from A to E
                (P(1) - lclLinkageNodeLoc(idxC, 1))^2 + ...
                (P(2) - lclLinkageNodeLoc(idxC, 2))^2 - y^2];   % dist from C to E
    options = optimoptions('fsolve', 'Display', 'none');  % Suppress output from fsolve
    PSol = fsolve(fun, initGuessE, options);
    lclLinkageNodeLoc(idxE, :) = [PSol(1), PSol(2), 0];
    initGuessE = [PSol(1), PSol(2)];
    theta(inc, 4) = atan2(lclLinkageNodeLoc(idxE, 2)-lclLinkageNodeLoc(originNode, 2), ...
                          lclLinkageNodeLoc(idxE, 1)-lclLinkageNodeLoc(originNode, 1));
    theta(inc, 5) = atan2(lclLinkageNodeLoc(idxE, 2)-lclLinkageNodeLoc(idxC, 2), ...
                          lclLinkageNodeLoc(idxE, 1)-lclLinkageNodeLoc(idxC, 1));

    if doFlip
        lclLinkageNodeLoc = mirrorY(lclLinkageNodeLoc, ...
                                    lclLinkageNodeLoc(originNode, :));
    end
    gblLinkageNodeLoc = local2globalcoord(lclLinkageNodeLoc', 'rr', lclOrigin', lclAxes);
    linkageNodeLoc(:, :, inc) = gblLinkageNodeLoc';
end

linkageNodeLoc(:, 3, :) = [];
individualDOFNodeLoc = linkageNodeLoc;
theta(theta < 0) = 2 * pi + theta(theta < 0);

if doFlip 
    dPhi = -rad2deg(diff(theta(:, [4,5])));
else
    dPhi = rad2deg(diff(theta(:, [4,5])));
end


end