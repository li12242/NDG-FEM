%======================================================================
%> @brief make uniform quadrilateral mesh object.
%>
%> @param N maximum polynomial order of basis functions
%> @param xlim limitation on x axis
%> @param ylim limitation on x axis
%> @param Mx Second argument
%> @param My Second argument
%> @param bcType Second argument
%>
%> @retval out1 return value for the first output variable
%> @retval out2 return value for the second output variable
%======================================================================
%> This function is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
%======================================================================
function [ mesh ] = makeUniformQuadMesh( N, xlim, ylim, Mx, My, bcType )
checkInput(xlim, ylim, Mx, My, bcType)

quad = StdQuad(N);
% Parameters
Nx = Mx + 1; % number of elements along x coordinate
Ny = My + 1;
K = Mx * My;
Nv = Nx * Ny;
EToR = double( enumRegion.Normal ) * ones(K, 1);

% Define vectex
% The vertex is sorted along x coordinate. (x coordinate counts first)

xmin = min(xlim); xmax = max(xlim);
ymin = min(ylim); ymax = max(ylim);
VX   = linspace(xmin, xmax, Nx) ;
VY   = linspace(ymin, ymax, Ny)';
VX   = repmat(VX, 1, Ny) ;
VY   = repmat(VY, 1, Nx)';
vx   = VX(:);
vy   = VY(:);

% Define EToV
% The element is conuting along x coordinate
EToV = zeros(4, Mx*My);
ind = 1;
for i = 1:My
    for j = 1:Mx
        % vertex index
        v1 = Nx*(i-1) + j;
        v2 = Nx*(i-1) + j + 1;
        v3 = Nx*i + j;
        v4 = Nx*i + j + 1;
        % Counterclockwise
        EToV(:, ind)=[v1, v2, v4, v3]';
        ind = ind +1;
    end% for
end% for

[ BCToV ] = makeUniformMeshBC( Mx, My, bcType );
mesh = NdgMesh2d( quad, Nv, vx, vy, K, EToV, EToR );
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
