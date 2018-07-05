function drawVerticalSlice( obj, cellId_2, nodeId_2, field3d )
%DRAWVERTICALSLICE Summary of this function goes here
%   Detailed explanation goes here

persistent p_handle

mesh3d = obj.mesh3d;
Nz = mesh3d.Nz;
Np2 = obj.mesh2d.cell.Np;
Npz = mesh3d.cell.Npz;
cellId_3 = (cellId_2 - 1) * Nz + (Nz:-1:1);
nodeId_3 = Np2 * ((1:Npz) - 1) + nodeId_2;

x = obj.mesh3d.x( nodeId_3, cellId_3 );
y = obj.mesh3d.y( nodeId_3, cellId_3 );
z = obj.mesh3d.z( nodeId_3, cellId_3 );

field = field3d( nodeId_3, cellId_3 );

plot3( x, y, z, 'k--' );
if ( isempty(p_handle) || ~isvalid(p_handle) )
    hold on;
    p_handle = plot3( x(:), y(:) + field(:), z(:), '-o' );
    box on;
    grid on;
else
    set( p_handle, 'YData', y(:) + field(:) );
end

end
