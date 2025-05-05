function [BC] = GetSupports(SUPP)
% Return degrees-of-freedom with fixed (prescribed) displacements
% Nf = sum(~isnan(SUPP(:, 2:3)), "all");
Nf = sum(SUPP(:, 2:3) == 1, "all");
BC = zeros(Nf, 1); j = 0;

for i = 1:size(SUPP, 1)
    if ~isnan(SUPP(i, 2)) && SUPP(i, 2) == 1, j = j + 1; BC(j) = 2 * SUPP(i) - 1; end
    if ~isnan(SUPP(i, 3)) && SUPP(i, 3) == 1, j = j + 1; BC(j) = 2 * SUPP(i);     end
end

if j~=Nf, error('Parsing number mismatch on BCs.'), end