function [ rhs ] = rhs_term(obj, f_Q )
%RHS_TERM_QUADRATURE Summary of this function goes here
%   Detailed explanation goes here

[ Eq, Gq ] = flux_term_quad(obj, f_Q); % volume flux term
[ num_fluxq ] = surf_term_quad( obj, f_Q ); % surface flux deviation

rhs = obj.mesh.cell.Drq'*(obj.mesh.rxq.*Eq) ...
    + obj.mesh.cell.Dsq'*(obj.mesh.sxq.*Eq) ...
    + obj.mesh.cell.Drq'*(obj.mesh.ryq.*Gq) ...
    + obj.mesh.cell.Dsq'*(obj.mesh.syq.*Gq) ...
    - obj.mesh.cell.Vbq'*( obj.mesh.fsq.*num_fluxq );

% for k = 1:obj.mesh.K
%     rhs(:,k) = obj.mesh.invM(:, :, k)*rhs(:,k);
% end

rhs = sum(bsxfun(@times, obj.mesh.invM, permute(rhs, [3,1,2])), 2);
rhs = permute(rhs, [1,3,2]);
end