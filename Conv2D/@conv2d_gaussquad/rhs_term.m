function [ rhs ] = rhs_term(obj, f_Q )
%RHS_TERM_QUADRATURE Summary of this function goes here
%   Detailed explanation goes here

[ Eq, Gq ] = flux_term_quad(obj, f_Q); % volume flux term
[ num_fluxq ] = surf_term_quad( obj, f_Q ); % surface flux deviation

rhs = zeros(obj.mesh.cell.Np, obj.mesh.K);
for k = 1:obj.mesh.K
rhs(:,k) = obj.mesh.Drq(:,:,k)*(obj.mesh.rxq(:, k).*Eq(:, k)) ...
    + obj.mesh.Dsq(:,:,k)*(obj.mesh.sxq(:, k).*Eq(:, k)) ...
    + obj.mesh.Drq(:,:,k)*(obj.mesh.ryq(:, k).*Gq(:, k)) ...
    + obj.mesh.Dsq(:,:,k)*(obj.mesh.syq(:, k).*Gq(:, k)) ...
    - obj.mesh.LIFTq(:,:,k)*( obj.mesh.fsq(:, k).*num_fluxq(:, k) );
end
end