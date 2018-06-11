function OutputOpenBoundaryVertCoor( obj )

filename = [ fileparts( mfilename('fullpath') ), ...
    '/tide/OpenBoundaryVertNodeCoordinate.txt' ];
fp = fopen( filename, 'w' );
for m1 = 1:obj.Nmesh
    
    edge = obj.meshUnion( m1 ).BoundaryEdge;
    edgeID = ( edge.ftype == enumBoundaryCondition.ClampedDepth );
    xb = edge.xb(:, edgeID);
    yb = edge.yb(:, edgeID);
    temp = [ xb(:)'; yb(:)' ];
    
    fprintf(fp, '%20.12f %20.12f\n', temp);
    
end

fclose(fp);

end
