%> @brief Initial function
%> This function will set the solver's mesh objects and initialize
%> the physical field
function initPhysFromOptions( obj, mesh )
% set mesh object
obj.meshUnion = mesh;
obj.Nmesh = numel(mesh);

% set edge object
[ obj.innerEdgeUnion, obj.haloEdgeUnion ] = setInnerEdgeUnion( obj, mesh );
obj.NhaloEdge = numel( obj.haloEdgeUnion );
% set option
obj.option = obj.setOption( obj.option );
% initilize physical field
obj.fphys = obj.setInitialField( );
end% func

function [ innerEdgeUnion, haloEdgeUnion ] = setInnerEdgeUnion( obj, mesh )
innerEdgeUnion = [];
for m = 1:obj.Nmesh
    if isa( mesh(m), 'NdgMesh2d' )
        edge = NdgInnerEdge2d( mesh(m), m );
        innerEdgeUnion = [ innerEdgeUnion, edge ];
    end
end

haloEdgeUnion = [];
for m1 = 1:obj.Nmesh
    for m2 = (m1+1):obj.Nmesh
        if isa( mesh(m1), 'NdgMesh2d' )
            edge = NdgHaloEdge2d( mesh(m1), mesh(m2), m1, m2 );
            haloEdgeUnion = [ haloEdgeUnion, edge ];
        end
    end
end

end
