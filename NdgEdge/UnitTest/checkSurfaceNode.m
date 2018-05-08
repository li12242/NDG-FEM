adv = AdvRotationHybridMesh2d( 1 );
mesh = adv.meshUnion;
edge1 = NdgHaloEdge2d( mesh, 1 );
edge2 = NdgHaloEdge2d( mesh, 2 );

% check edge1 surface nodes
edge = edge2;
xM = zeros(edge.Nfp, edge.Ne);
yM = zeros(edge.Nfp, edge.Ne);
xP = zeros(edge.Nfp, edge.Ne);
yP = zeros(edge.Nfp, edge.Ne);

for i = 1:edge.Ne
    m1 = edge.FToM(1, i);
    m2 = edge.FToM(2, i);
    
    xM(:, i) = mesh(m1).x(edge.FToN1(:, i), edge.FToE(1, i));
    yM(:, i) = mesh(m1).y(edge.FToN1(:, i), edge.FToE(1, i));
    
    xP(:, i) = mesh(m2).x(edge.FToN2(:, i), edge.FToE(2, i));
    yP(:, i) = mesh(m2).y(edge.FToN2(:, i), edge.FToE(2, i));
end

xM - xP
yM - yP


