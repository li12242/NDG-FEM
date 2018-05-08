function [ meshUnion ] = makeMeshUnion( Nmesh, varargin )

if nargin ~= (Nmesh + 1)
    msgID = [ mfilename, ':InputMeshNumberError'];
    msgtext = 'The number of input meshes is not equal to Nmesh.';
    throw( MException(msgID, msgtext) );
end

meshUnion = [];
for m = 1:Nmesh
    meshUnion = [ meshUnion, varargin{m} ];
end

for m = 1:Nmesh
    mesh = meshUnion(m);
    for m1 = 1:Nmesh
        if( m == m1 ) continue; end
        mesh1 = meshUnion(m1);
        mesh.assembleMeshConnection(mesh1, m, m1);
    end
end

% for m = 1:Nmesh
%     mesh = meshUnion(m);
%     for m1 = 1:Nmesh
%         if( m == m1 ) continue; end
%         mesh1 = meshUnion(m1);
%         mesh.assembleNdgEdgeConnection(mesh1, m, m1);
%     end
% end

end

