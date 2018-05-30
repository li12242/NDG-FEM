function FindOBVertex( obj )

for m1 = 1:obj.Nmesh
    mesh = obj.meshUnion( m1 );
    b = mesh.K;
    a = StdQuad.Nv;
    
    edgeID = (mesh.EToB == NdgEdgeType.ClampedDepth);
    [f, k] = find(edgeID);
    vid = mesh.cell.FToV(:, f);
    kid = [k, k]';
    ind = sub2ind([a,b], vid(:), kid(:));
    
    obj.OBVid(:,1) = unique( mesh.EToV(ind) );
    obj.OBVid(:,2) = mesh.vx(obj.OBVid(:,1));
    obj.OBVid(:,3) = mesh.vy(obj.OBVid(:,1));
    
end

end

