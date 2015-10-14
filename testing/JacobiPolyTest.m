function error = JacobiPolyTest(N)
%% JacobiP
nModePoints = N + 1;
[nodalPoints,~] = Polylib.zwglj(nModePoints);
V1D = zeros(length(nodalPoints),nModePoints);
for j=1:nModePoints
    V1D(:,j) = JacobiP(nodalPoints(:), 0, 0, j-1);
end
%% Polylib
for j=1:nModePoints
    V2D(:,j) = Polylib.jacobfd(nodalPoints, j-1);
end
%% compare
error = abs(V1D - V2D);
end