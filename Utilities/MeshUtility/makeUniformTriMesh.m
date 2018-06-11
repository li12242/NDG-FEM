function [ mesh ] = makeUniformTriMesh( N, xlim, ylim, Mx, My, bcType )
checkInput(xlim, ylim, Mx, My, bcType);
tri = StdTri(N);

flag = 0;
% Parameters
Nx = Mx + 1; % number of nodes along x coordinate
Ny = My + 1;
K  = Mx * My * 2;
Nv = Nx * Ny;
EToR = enumRegion.Normal * ones(K, 1, 'int8');

% Define vertex
% The vertex is sorted along x coordinate. (x coordinate counts first)
xmin = min( xlim ); xmax = max( xlim );
ymin = min( ylim ); ymax = max( ylim );
VX   = linspace(xmin, xmax, Nx) ;
VY   = linspace(ymin, ymax, Ny)'; 
VX   = repmat(VX, 1, Ny) ;
VY   = repmat(VY, 1, Nx)'; 
vx   = VX(:);
vy   = VY(:);

% Define EToV
% The element is conuting along x coordinate
EToV = zeros(3, 2*Mx*My);
for i = 1:My % each row
    for j = 1:Mx
        % element index
        ind1 = 2*Mx*(i-1) + j;
        ind2 = 2*Mx*(i-1)+Mx+j;
        % vertex index
        v1 = Nx*(i-1) + j;
        v2 = Nx*(i-1) + j + 1;
        v3 = Nx*i + j;
        v4 = Nx*i + j + 1;
        % Counterclockwise
        if flag % '/' divided
            EToV(:, ind1) = [v1, v4, v3]';
            EToV(:, ind2) = [v1, v2, v4]';
        else    % '\' divided
            EToV(:, ind1) = [v1, v2, v3]';
            EToV(:, ind2) = [v2, v4, v3]';
        end% if
    end
end

[ BCToV ] = makeUniformMeshBC( Mx, My, bcType );
mesh = NdgMesh2d( tri, Nv, vx, vy, K, EToV, EToR );
mesh.ConnectMeshUnion( 1, mesh);
mesh.InnerEdge = NdgInnerEdge2d( mesh, mesh.ind );
mesh.BoundaryEdge = NdgHaloEdge2d( mesh, mesh.ind, BCToV );
end

function checkInput(xlim, ylim, Mx, My, bcType)
if (numel(xlim) ~= 2) || (numel(ylim) ~= 2)
    msgID = 'makeUniformTriUMeshUnion:InputLimsError';
    msgtext = 'The input xlim and ylim should be a veoctor with 2 values.';
    ME = MException(msgID, msgtext);
    throw(ME);
end

if (Mx <=0 ) || (My <= 0 )
    msgID = 'makeUniformTriUMeshUnion:InputCellNumError';
    msgtext = 'The input Mx and My should be a integer greater than 0.';
    ME = MException(msgID, msgtext);
    throw(ME);
end

if (numel(bcType) ~= 4)
    msgID = 'makeUniformTriUMeshUnion:InputBoundaryConditionError';
    msgtext = 'The input bcType should contain four boundary types.';
    ME = MException(msgID, msgtext);
    throw(ME);
end

end

