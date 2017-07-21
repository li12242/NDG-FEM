function [ Eq, Gq ] = flux_term_quad( obj, f_Q )
%FLUX_TERM_QUAD Summary of this function goes here
%   Detailed explanation goes here

fq_Q = obj.mesh.cell.project_node2quad(f_Q);
Eq = obj.uq.*fq_Q;
Gq = obj.vq.*fq_Q;
end

