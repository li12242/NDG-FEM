function tests = TestSlopeLimitBiswas
tests = functiontests(localfunctions);
end% func

function testIncreasing(testCase)
N = 2;
line = StdRegions.Line(N); EToV = [1,2; 2,3; 3,4]; VX = [1,2,3,4];
mesh = MultiRegions.RegionLine(line, EToV, VX);

u = mesh.x;
ulimit = Utilities.Limiter.SlopeLimitBiswas(mesh, u);

expSolution = u(:,2:3);
verifyEqual(testCase, ulimit(:,2:3), expSolution)
end% func

% function testDecreasing
% end% func