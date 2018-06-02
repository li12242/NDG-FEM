function [ Nfp, M ] = assembleMassMatrix( obj )
    cell = StdLine( obj.mesh.cell.N );
    M = cell.M;
    Nfp = cell.Np;
end