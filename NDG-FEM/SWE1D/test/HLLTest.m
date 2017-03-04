function tests = HLLTest
tests = functiontests(localfunctions);
end% func

function testLeftIsDry1(testCase)
% hs > 0
hM = 0; hP = 4; uM = 0; uP = 2; mesh.nx = 1; g = 9.8;
[SM, SP] = EstimateWaveSpeed(mesh, hM, hP, uM, uP);
expSolution = [uP - 2*sqrt(g*hP), uP + sqrt(g*hP)];
verifyEqual(testCase,[SM, SP],expSolution)
end% func

function testLeftIsDry2(testCase)
% hs = 0
hM = 0; hP = 4; uM = 0; uP = 14; mesh.nx = 1; g = 9.8;
[SM, SP] = EstimateWaveSpeed(mesh, hM, hP, uM, uP);
expSolution = [uP - 2*sqrt(g*hP), uP + sqrt(g*hP)];
verifyEqual(testCase,[SM, SP],expSolution)
end% func


function testRightIsDry1(testCase)
% hs > 0
hM = 10; hP = 0; uM = 1; uP = 0; mesh.nx = 1; g = 9.8;
[SM, SP] = EstimateWaveSpeed(mesh, hM, hP, uM, uP);
expSolution = [uM - sqrt(g*hM), uM + 2*sqrt(g*hM)];
verifyEqual(testCase,[SM, SP],expSolution)
end% func

function testRightIsDry2(testCase)
% hs = 0
hM = 10; hP = 0; uM = -20; uP = 0; mesh.nx = 1; g = 9.8;
[SM, SP] = EstimateWaveSpeed(mesh, hM, hP, uM, uP);
expSolution = [uM - sqrt(g*hM), uM + 2*sqrt(g*hM)];
verifyEqual(testCase,[SM, SP],expSolution)
end% func

function testShockWave(testCase)

end% func

function testRafactionWave(testCase)

end% func