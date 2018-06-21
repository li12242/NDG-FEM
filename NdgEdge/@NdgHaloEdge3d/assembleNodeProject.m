function [ obj ] = assembleNodeProject( obj, mesh )
%ASSEMBLENODEPROJECT Summary of this function goes here
%   Detailed explanation goes here

mesh1 = obj.mesh;
Nfp = obj.Nfp;
FToN1 = zeros( Nfp, obj.Ne );
FToN2 = zeros( Nfp, obj.Ne );
nx = zeros( Nfp, obj.Ne );
ny = zeros( Nfp, obj.Ne );
nz = zeros( Nfp, obj.Ne );
Js = zeros( Nfp, obj.Ne );
Fmask1 = mesh1.cell.Fmask;
Nfp = mesh1.cell.Nfp;

for n = 1:obj.Ne
    m2 = obj.FToM(2, n);
    e1 = obj.FToE(1, n);
    e2 = obj.FToE(2, n);
    f1 = obj.FToF(1, n);
    f2 = obj.FToF(2, n);
    Fmask2 = mesh(m2).cell.Fmask;

    % set local node index
    FToN1(:, n) = Fmask1(:, f1);

    x2 = mesh(m2).x(Fmask2(:, f2), e2);
    y2 = mesh(m2).y(Fmask2(:, f2), e2);
    z2 = mesh(m2).z(Fmask2(:, f2), e2);
    for m = 1 : Nfp
        ind = FToN1(m, n);
        x1 = mesh1.x(ind, e1);
        y1 = mesh1.y(ind, e1);
        z1 = mesh1.z(ind, e1);

        t = ( abs( x2 - x1 ) + abs( y2 - y1 ) + abs( z2 - z1 ) ) < 1e-10;
        FToN2(m, n) = Fmask2(t, f2);
    end

    % set outward normal vector
    if mesh1.cell.type == enumStdCell.PrismTri
        [ nx( :, n ), ny( :, n ), nz( :, n ), Js( :, n ) ] ...
            = obj.PrismTriJacobian3d( mesh1, f1, e1, Fmask1( Nfp(f1), f1 ));
    elseif mesh1.cell.type == enumStdCell.PrismQuad
        [ nx( :, n ), ny( :, n ), nz( :, n ), Js( :, n ) ] ...
            = obj.PrismQuadJacobian3d( mesh1, f1, e1, Fmask1( Nfp(f1), f1 ));
    end
end

obj.FToN1 = FToN1;
obj.FToN2 = FToN2;
obj.nx = nx;
obj.ny = ny;
obj.nz = nz;
obj.Js = Js;

ind = obj.FToN1 + mesh1.cell.Np * repmat( obj.FToE(1, :) - 1, obj.Nfp, 1 );
obj.xb = mesh1.x( ind );
obj.yb = mesh1.y( ind );
obj.zb = mesh1.z( ind );

end

