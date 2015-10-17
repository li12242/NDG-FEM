function tests = TestSlopeLimitSWE
tests = functiontests(localfunctions);
end% func

function testIncreasing(testCase)
[EToV, VX, VY, EToR, BC] = ...
    Utilities.Mesh.MeshReaderTriangle('testing/Utilities/Limiter/mesh/simple');

N = 2;
tri = StdRegions.Triangle(N);
mesh = MultiRegions.RegionTri(tri, EToV, VX, VY);

h = mesh.x;
hlimit = Utilities.Limiter.SlopeLimitSWE(mesh, h);

verifyEqual(testCase, hlimit, mesh.x)
end% func

function testShock(testCase)
[EToV, VX, VY, EToR, BC] = ...
    Utilities.Mesh.MeshReaderTriangle('testing/Utilities/Limiter/mesh/shock');
N = 2;
tri = StdRegions.Triangle(N);
mesh = MultiRegions.RegionTri(tri, EToV, VX, VY);

xc = mean(mesh.x); left = xc < 0.5; right = xc > 0.5;
h = zeros(size(mesh.x));
h(:, left) = 0; h(:, right) = 1;

hlimit = Utilities.Limiter.SlopeLimitSWE(mesh, h);

verifyEqual(testCase, hlimit, mesh.x)
end% func