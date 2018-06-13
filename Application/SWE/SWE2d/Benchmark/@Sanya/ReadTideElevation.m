function ReadTideElevation( obj )
%READTIDEELEVATION
data = load( obj.tidalFile );

TotalNb = 0;
Nb = cell( obj.Nmesh, 1 );
obj.OBEdgeIndex = cell( obj.Nmesh, 1 );
for m = 1:obj.Nmesh
    edge = obj.meshUnion( m ).BoundaryEdge;
    obj.OBEdgeIndex{m} = find( edge.ftype == enumBoundaryCondition.ClampedDepth );
    Nb{m} = edge.Nfp * numel( obj.OBEdgeIndex{m} );
    TotalNb = TotalNb + Nb{m};
end
Ntime = numel( data ) / TotalNb;

obj.Tide = cell( obj.Nmesh, 1 );
for m = 1:obj.Nmesh
    Nbedge = numel( obj.OBEdgeIndex{m} );
    obj.Tide{m} = zeros( edge.Nfp, Nbedge, Ntime );
end

counter = 1;
for t = 1:Ntime
    for m = 1:obj.Nmesh
        mesh = obj.meshUnion( m );
        edge = obj.meshUnion( m ).BoundaryEdge;
        
        ind = edge.FToN1 + repmat(edge.FToE(1, :) - 1, edge.Nfp, 1) * mesh.cell.Np;
        bot = obj.fphys{m}(:, :, 4);
        obj.fext{m}(:, :, 4) = bot( ind );
        
        Nbedge = numel( obj.OBEdgeIndex{m} );
        obj.Tide{m}(:, :, t) = reshape( ...
            data(counter:(counter + Nb{m}-1)), edge.Nfp, Nbedge );
        
        counter = counter + Nb{m};
    end
end

end
