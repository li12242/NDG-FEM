%======================================================================
%> @brief make the boundary conditions for the 2d uniform mesh
%>
%> More detailed description.
%>
%> @param Mx number of elements on x axis;
%> @param My number of elements on y axis;
%> @param face_type boundary types for north, sourth, west, east boundaries;
%>
%> @retval BCToV boundary types and the index of vertices on edges
%======================================================================
%> This function is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
%======================================================================
function [ BCToV ] = makeUniformMeshBC( Mx, My, bcType )

BCToV = ones((Mx+My)*2, 3);
Nx = Mx + 1; % nodes on x axis
Ny = My + 1; % nodes on y axis
% sourth boundary
st = 1; 
se = st + Mx - 1;
BCToV(st:se,[1,2]) = [1:(Nx-1); 2:Nx]';
BCToV(st:se, 3) = bcType(1);
% north boundary
st = se + 1; 
se = st + Mx - 1;
vs = Nx*(Ny-1)+1; 
ve = Nx*Ny;
BCToV(st:se,[1,2]) = [vs:(ve-1); (vs+1):ve]';
BCToV(st:se, 3) = bcType(2);
% west boundary
st = se + 1; 
se = st + My - 1;
vs = 1; 
ve = Nx*(Ny-1)+1;
vstrid = Nx;
BCToV(st:se,[1,2]) = [vs:vstrid:(ve-vstrid); (vs+vstrid):vstrid:ve]';
BCToV(st:se, 3) = bcType(3);
% east boundary
st = se + 1; 
se = st + My - 1;
vs = Nx; 
ve = Nx*Ny; 
vstrid = Nx;
BCToV(st:se,[1,2]) = [vs:vstrid:(ve-vstrid); (vs+vstrid):vstrid:ve]';
BCToV(st:se, 3) = bcType(4);

BCToV = BCToV';
end

