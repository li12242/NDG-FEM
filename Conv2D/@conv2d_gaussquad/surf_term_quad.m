function [ num_flux ] = surf_term_quad( obj, f_Q )
%SURF_TERM_QUAD Summary of this function goes here
%   Detailed explanation goes here

ubq = obj.mesh.cell.project_node2surf_quad( obj.u );
vbq = obj.mesh.cell.project_node2surf_quad( obj.v );
f_extbq = obj.mesh.cell.project_node2surf_quad( obj.f_extQ );
f_Qbq = obj.mesh.cell.project_node2surf_quad( f_Q );

num_flux = upwind_flux(f_Qbq, f_extbq, ubq, vbq, ...
    obj.mesh.nxq, obj.mesh.nyq, obj.mesh.eidMq, obj.mesh.eidPq, ...
    obj.mesh.eidtypeq);
end

