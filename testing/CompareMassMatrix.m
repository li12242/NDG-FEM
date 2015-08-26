function error = CompareMassMatrix(N)
%% Gauss integrate

nModePoints = N+1;
nQuadPoints = N+3;
[nodalPoints,~] = Polylib.zwglj(nModePoints);
[quadPoints, quadWeights] = Polylib.zwglj(nQuadPoints);
base = zeros(nModePoints, nQuadPoints);
for i= 1:nModePoints
    for j = 1:nQuadPoints
        base(i,j) = Polylib.hglj(i-1,quadPoints(j),nodalPoints,nModePoints);
    end
end

weight = zeros(nQuadPoints,nQuadPoints);
for j = 1:nQuadPoints
    weight(j,j) = quadWeights(j);
end
M_1 = base*weight*base';

%% Polylib Vandmonde
for j=1:nModePoints
    [p, ~] = Polylib.jacobfd(nodalPoints, j-1);
    V2D(:,j) = p(:)*sqrt((2*(j-1)+1)/2);
end
M_3 = inv(V2D*V2D');

%% Vandmode Matrix
% Compute basic Legendre Gauss Lobatto grid
r = JacobiGL(0,0,N);

% Build reference element matrices
V  = Vandermonde1D(N, r);

M_2 = inv(V*V');

%% The 

error = abs(M_1 - M_3);
end