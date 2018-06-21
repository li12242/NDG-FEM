function obj = assembleNodeProject( obj, mesh )

% reference element
cell = mesh.cell;

Nfp = obj.Nfp;
Ne = obj.Ne;
FToN1 = zeros( Nfp, Ne );
FToN2 = zeros( Nfp, Ne );
nx = zeros( Nfp, Ne ); ny = zeros( Nfp, Ne );
nz = zeros( Nfp, Ne ); Js = zeros( Nfp, Ne );
Fmask = cell.Fmask;
% [ r, s ] = assembleInpterpNode( Nh, Nz );
% [ Dr, Ds ] = nodal_derivative_func( Nh, Nz, r, s );

for n = 1 : Ne
    e1 = obj.FToE(1, n);
    e2 = obj.FToE(2, n);
    f1 = obj.FToF(1, n);
    f2 = obj.FToF(2, n);

    % set local node index
    FToN1(:, n) = Fmask(:, f1);

    x2 = mesh.x(Fmask(:, f2), e2);
    y2 = mesh.y(Fmask(:, f2), e2);
    z2 = mesh.z(Fmask(:, f2), e2);
    for m = 1 : Nfp
        ind = FToN1(m, n);
        x1 = mesh.x(ind, e1);
        y1 = mesh.y(ind, e1);
        z1 = mesh.z(ind, e1);

        t = ( abs( x2 - x1 ) + abs( y2 - y1 ) + abs( z2 - z1 ) ) < 1e-10;
        FToN2(m, n) = Fmask(t, f2);
    end

    if mesh.cell.type == enumStdCell.PrismTri
        [ nx( :, n ), ny( :, n ), nz( :, n ), Js( :, n ) ] ...
            = obj.PrismTriJacobian3d( mesh, f1, e1, cell.Fmask(1:cell.Nfp(f1), f1));
    elseif mesh.cell.type == enumStdCell.PrismQuad
        [ nx( :, n ), ny( :, n ), nz( :, n ), Js( :, n ) ] ...
            = obj.PrismQuadJacobian3d( mesh, f1, e1, cell.Fmask(1:cell.Nfp(f1), f1));
    end
end

obj.FToN1 = FToN1;
obj.FToN2 = FToN2;
obj.nx = nx; 
obj.ny = ny; 
obj.nz = nz;
obj.Js = Js;

end% func
