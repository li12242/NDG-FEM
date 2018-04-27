function [ FToN1, FToN2, nx, ny, nz, Js ] = assembleNodeProject( obj, mesh )

cell = mesh.cell;
FToN1 = zeros( obj.bcell.Np, obj.Ne );
FToN2 = zeros( obj.bcell.Np, obj.Ne );
nx = zeros( obj.bcell.Np, obj.Ne );
ny = zeros( obj.bcell.Np, obj.Ne );
nz = zeros( obj.bcell.Np, obj.Ne );
Js = zeros( obj.bcell.Np, obj.Ne );

for n = 1:obj.Ne
    k1 = obj.FToE(1, n);
    k2 = obj.FToE(2, n);
    f1 = obj.FToF(1, n);
    f2 = obj.FToF(2, n);
    vert1 = mesh.EToV( cell.FToV(:, f1), k1);
    vert2 = mesh.EToV( cell.FToV(:, f2), k2);
    
    % set local node index
    FToN1(:, n) = mesh.cell.Fmask(:, f1);
    if vert2(1) == vert1(1)
        FToN2(:, n) = mesh.cell.Fmask(:, f2);
    else
        FToN2(:, n) = flip( mesh.cell.Fmask(:, f2) );
    end
    
    % set outward normal vector
    tmp = sum( cell.Nfp(1:f1) );
    fid = (tmp - cell.Nfp(f1) + 1):tmp;
    nx(:, n) = mesh.nx(fid, k1);
    ny(:, n) = mesh.ny(fid, k1);
    nz(:, n) = mesh.nz(fid, k1);
    Js(:, n) = mesh.Js(fid, k1);
end

end

