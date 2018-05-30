%> Update the external field
function matUpdateExternalField( obj, time, fphys )
matUpdateExternalField@NdgPhysMat( obj, time, fphys );

for m1 = 1:obj.Nmesh
    mesh = obj.meshUnion( m1 );
    delta = obj.tideinterval;   
    
    if abs(time-0) < 1e-6
        tidet = obj.Tide(:,1);
    else
        a = ceil(time/delta); 
        tide1 = obj.Tide(:,a);
        tide2 = obj.Tide(:,a+1);
        tidet = (time/delta-a+1)*(tide2-tide1)+tide1;
    end
    
    fvT = zeros(mesh.Nv,1);
    for i = 1 : numel(obj.OBVid(:,1))
        fvT(obj.OBVid(i,1)) = tidet(i);
    end
    ele_fvT = fvT(mesh.EToV);
    fnT = mesh.cell.project_vert2node(ele_fvT);
    
    ele = fnT - obj.fphys{m1}(:,:,4);
    obj.fext{m1}(:,:,1) = max(0,ele);
end

% obj.FextLoader.getExternalField(obj,time);

end% func