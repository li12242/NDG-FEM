function testing_Line
%% testing line basic
nOrder = 8;
point = StdRegions.Point;
line = StdRegions.LineBasic(nOrder);

% Compare Mass Matrix
TOL = 10e-10;
nModePoints = nOrder+1;
[nodalPoints,~] = Polylib.zwglj(nModePoints);
for j=1:nModePoints
    [p, ~] = Polylib.jacobfd(nodalPoints, 0,0, j-1);
    V2D(:,j) = p(:)*sqrt((2*(j-1)+1)/2);
end
M_3 = inv(V2D*V2D');

if((line.M - M_3)< TOL) 
    fprintf('Mass Matrix Correct\n');
end
% Compare Derivative Matrix
D = Polylib.Dglj(nodalPoints); 
D = D';

if((line.Dr - D)< TOL) 
    fprintf('Derivative Matrix Correct\n');
end

%% testing line element
nOrder = 8;
line = StdRegions.Line(nOrder);
% testing face node list
nodelist = line.getFaceListToNodeList; y = zeros(size(line.r(nodelist)));
figure; hold on;
plot(line.r, line.r.*0, 'ro');
plot(line.r(nodelist), y, 'bo')
for i = 1:line.nFaceNode
    text(line.r(nodelist(i))+0.05, y(i)+0.05, num2str(nodelist(i)))
end
% testing edge mass matrix


end