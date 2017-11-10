function [ rhs ] = rhs_term(obj, f_Q )
%RHS_TERM Calculate the right hand term
%   Detailed explanation goes here

[ Eq, Gq ] = flux_term_quad(obj, f_Q); % volume flux term
[ num_fluxq ] = surf_term_quad( obj, f_Q ); % surface flux deviation

rhs = obj.mesh.cell.Drq'*(obj.mesh.rxwJq.*Eq) ...
    + obj.mesh.cell.Dsq'*(obj.mesh.sxwJq.*Eq) ...
    + obj.mesh.cell.Drq'*(obj.mesh.rywJq.*Gq) ...
    + obj.mesh.cell.Dsq'*(obj.mesh.sywJq.*Gq) ...
    - obj.mesh.cell.Vbq'*( obj.mesh.wJsq.*num_fluxq );

% for k = 1:obj.mesh.K
%     rhs(:,k) = obj.mesh.invM(:, :, k)*rhs(:,k);
% end

rhs = sum(bsxfun(@times, obj.mesh.invM, permute(rhs, [3,1,2])), 2);
rhs = permute(rhs, [1,3,2]);
end