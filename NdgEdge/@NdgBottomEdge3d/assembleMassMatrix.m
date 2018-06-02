function [ Nfp, M ] = assembleMassMatrix( obj, N, Nz )
    
    if ( obj.mesh.cell.type == NdgCellType.PrismTri )
        cell = StdTri(N);
        Nfp = cell.Np;
        M = cell.M;
    elseif( obj.mesh.cell.type == NdgCellType.PrismQuad )
        cell = StdQuad(N);
        Nfp = cell.Np;
        M = cell.M;
    end
end