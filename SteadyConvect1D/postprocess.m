ne = [10, 20, 40, 80];

% get mesh
x1 = 0; x2 = 2*pi; % domain

N = 3; %nElement = 20;
BC = [2,1];
line = StdRegions.Line(N);

L1 = zeros(size(ne)); L2 = zeros(size(ne));
for i = 1:numel(ne)
    [~, VX, ~, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, ne(i));
    mesh = MultiRegions.RegionLineBC(line, EToV, VX, BC);
    
    
    u = SteadyConvectionDriver(N, ne(i));
    exu = -cos(mesh.x) + 2;
    L1(i) = sum(sum(abs(u - exu)))./mesh.nNode;
    L2(i) = sqrt( sum( sum( (u - exu).^2./mesh.nNode ) ) );
end% for

rate1 = zeros(numel(ne) - 1, 1);
rate2 = zeros(numel(ne) - 1, 1);

fprintf('| ne | L1 | rate |\n')
fprintf('| --- | --- | --- |\n')
fprintf('| %d | %f | \\ |\n', ne(1), L1(1));
for i = 1:numel(ne)-1
    rate1(i) = log2( L1(i)./L1(i+1) );
    fprintf('| %d | %f | %f |\n', ne(i+1), L1(i+1), rate1(i));
end% for

fprintf('| ne | L2 | rate |\n')
fprintf('| --- | --- | --- |\n')
fprintf('| %d | %f | \\ |\n', ne(1), L2(1));
for i = 1:numel(ne)-1
    rate2(i) = log2( L2(i)./L2(i+1) );
    fprintf('| %d | %f | %f |\n', ne(i+1), L2(i+1), rate2(i));
end% for

