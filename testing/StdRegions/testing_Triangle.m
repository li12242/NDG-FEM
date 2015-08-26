function testing_Triangle
%% testing triangle basic element
TOTALERR = 10^-10;
nOrder = 5;
tri = StdRegions.TriangleBasic(nOrder);
% testing mass matrix
[x,y] = Nodes2D(nOrder); [r,s] = xytors(x,y);
V = Vandermonde2D(nOrder,r,s); invV = inv(V);
MassMatrix = invV'*invV;
if sum(abs(tri.M - MassMatrix) > TOTALERR)
    fprintf('Mass Matrix error occurs\n')
else
    fprintf('Mass Matrix \t testing pass\n')
end

% testing deri matrix
[Dr,Ds] = Dmatrices2D(nOrder, r, s, V);
if sum(abs(tri.Dr - Dr) > TOTALERR)
    fprintf('Derivative Matrix of r error occurs\n')
else
    fprintf('Derivative Matrix of r \t testing pass\n')
end

if sum(abs(tri.Ds - Ds) > TOTALERR)
    fprintf('Derivative Matrix of s error occurs\n')
else
    fprintf('Derivative Matrix of s \t testing pass\n')
end

fprintf('LineBasic \t testing past\n')
fprintf('\n')

% fprintf('press Enter to continue\n')
% pause

%% tesing triangle element
nOrder = 5;
tri = StdRegions.Triangle(nOrder);

% test edge node list
nodelist = tri.getFaceListToNodeList;
figure; hold on;
plot(tri.r, tri.s, 'ro');
plot(tri.r(nodelist), tri.s(nodelist), 'bo')
for i = 1:tri.nFaceNode
    text(tri.r(nodelist(i))+0.05, tri.s(nodelist(i))+0.05, num2str(nodelist(i)))
end

vorder = [2,1];
for iface = 1:3
    facelist = tri.getReorderFaceListAtFace(iface, vorder);
    nodelist(facelist)
end
% test edge Mass Matrix
% the edge mass matrix is made of mass matrix of face element -line
for iface = 1:tri.nFace
    [~,nodelist] = tri.getNodeListAtFace(iface);
    facelist = tri.getFaceListAtFace(iface);
    temp = tri.FaceMassMatrixSmall(nodelist, facelist);
    if any(temp - line.M)
        error('properties: FaceMassMatrixSmall error occurs');
    end
end
% check full edge mass matrix
vmapM = tri.getFaceListToNodeList;
tri.FaceMassMatrixFull*tri.r - tri.FaceMassMatrixSmall*tri.r(vmapM)

end