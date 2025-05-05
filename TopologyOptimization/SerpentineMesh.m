function [NODE, ELEM, SUPP, LOAD] = SerpentineMesh(Nx, Ny, Lx, Ly)
NODE = zeros((Nx + 1) * (Ny + 1), 2);
X = linspace(0, Lx, Nx + 1);
syms u v
circ1(u, v) = (u - 2)^2 + (v - 9.9)^2 - (10.1)^2;
circ2(u, v) = (u - 2)^2 + (v - 10.9)^2 - (10.1)^2;
circ3(u, v) = (u - 6)^2 + (v + 9.9)^2 - (10.1)^2;
circ4(u, v) = (u - 6)^2 + (v + 8.9)^2 - (10.1)^2;
range = [-(Ly+0.5), Ly+0.5];
for i = 1:length(X)
    if X(i) <= Lx/2
        NODE(2*i - 1, :) = [X(i), vpasolve(circ1(X(i), v), v, range)];
        NODE(2*i, :) = [X(i), vpasolve(circ2(X(i), v), v, range)];
    else
        NODE(2*i - 1, :) = [X(i), vpasolve(circ3(X(i), v), v, range)];
        NODE(2*i, :) = [X(i), vpasolve(circ4(X(i), v), v, range)];
    end
end

k = 0;
for j=1:Ny, for i=1:Nx
        k = k+1;
        n1 = (i-1)*(Ny+1)+j; n2 = i*(Ny+1)+j;
        ELEM{k} = [n1 n2 n2+1 n1+1];
end, end

SUPP = [1, 1, 0;
        Ny + 1, 1, 1];
LOAD = [Nx * (Ny + 1) + 1, 0, -1];
end