function obj = assembleNodeProject( obj, mesh )

mesh = obj.mesh;
cell = mesh.cell;
Nfp = obj.Nfp;
Ne = obj.Ne;
FToN1 = zeros( Nfp, Ne );
FToN2 = zeros( Nfp, Ne );
nx = zeros( Nfp, Ne );
ny = zeros( Nfp, Ne );
nz = zeros( Nfp, Ne );
Js = zeros( Nfp, Ne );

for n = 1:obj.Ne
    k1 = obj.FToE(1, n);
    k2 = obj.FToE(2, n);
    f1 = obj.FToF(1, n);
    f2 = obj.FToF(2, n);
    vert1 = mesh.EToV( cell.FToV(:, f1), k1);
    vert2 = mesh.EToV( cell.FToV(:, f2), k2);
    
    % set local node index
    FToN1(:, n) = cell.Fmask(:, f1);
    if vert2(1) == vert1(1)
        FToN2(:, n) = mesh.cell.Fmask(:, f2);
    else
        FToN2(:, n) = flip( mesh.cell.Fmask(:, f2) );
    end
    
    if mesh.cell.type == enumStdCell.Tri
        [ nx(:, n), ny(:, n), Js(:, n) ] = TriJacobian2d( mesh, f1, k1, FToN1(:, n) );
    elseif mesh.cell.type == enumStdCell.Quad
        [ nx(:, n), ny(:, n), Js(:, n) ] = QuadJacobian2d( mesh, f1, k1, FToN1(:, n) );
    end

end

obj.FToN1 = FToN1;
obj.FToN2 = FToN2;
obj.nx = nx;
obj.ny = ny;
obj.nz = nz;
obj.Js = Js;

end

function [ nx, ny, Js ] = QuadJacobian2d( mesh, f1, e1, nodeId )
rx = mesh.rx( nodeId, e1 ); ry = mesh.ry( nodeId, e1 );
sx = mesh.sx( nodeId, e1 ); sy = mesh.sy( nodeId, e1 );

if f1 == 1
    nx = - sx;
    ny = - sy;
elseif f1 == 2
    nx = rx;
    ny = ry;
elseif f1 == 3
    nx = sx;
    ny = sy;
elseif f1 == 4
    nx = - rx;
    ny = - ry;
end

Js = sqrt( nx .* nx + ny .* ny );
nx = nx ./ Js;
ny = ny ./ Js;
Js = Js .* mesh.J( nodeId, e1 );

end

function [ nx, ny, Js ] = TriJacobian2d( mesh, f1, e1, nodeId )
rx = mesh.rx( nodeId, e1 ); ry = mesh.ry( nodeId, e1 );
sx = mesh.sx( nodeId, e1 ); sy = mesh.sy( nodeId, e1 );

if f1 == 1
    nx = - sx;
    ny = - sy;
elseif f1 == 2
    nx = rx + sx;
    ny = ry + sy;
elseif f1 == 3
    nx = - rx;
    ny = - ry;
end

Js = sqrt( nx .* nx + ny .* ny );
nx = nx ./ Js;
ny = ny ./ Js;
Js = Js .* mesh.J( nodeId, e1 );

end

