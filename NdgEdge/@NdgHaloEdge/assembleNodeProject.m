function [ IntM, FToN1, FToN2 ] = assembleNodeProject( obj, mesh1, mesh2 )

bcell1 = obj.bcell1;
bcell2 = obj.bcell2;
IntM = zeros( bcell1.Np, bcell2.Np );
for n = 1:bcell2.Np
    IntM(:, n) = bcell2.orthogonal_func( bcell2.N, n, bcell1.r, bcell1.s, 0);
end
IntM = IntM/bcell2.V;

cell1 = mesh1.cell;
cell2 = mesh2.cell;
FToN1 = zeros( bcell1.Np, obj.M );
FToN2 = zeros( bcell2.Np, obj.M );

for n = 1:obj.M
    k1 = obj.FToE(1, n);
    k2 = obj.FToE(2, n);
    f1 = obj.FToF(1, n);
    f2 = obj.FToF(2, n);
    vert1 = mesh1.EToV( cell1.FToV(:, f1), k1);
    vert2 = mesh2.EToV( cell2.FToV(:, f2), k2);
    
    % set local node index
    FToN1(:, n) = mesh1.cell.Fmask(:, f1);
    % set adjacent node index
    if vert2(1) == vert1(1)
        FToN2(:, n) = mesh2.cell.Fmask(:, f2);
    else
        FToN2(:, n) = flip( mesh2.cell.Fmask(:, f2) );
    end
end

end% func