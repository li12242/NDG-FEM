function cP = nei_node_val(obj, cM, cP, c_extM, ftype)
%ADJ_NODE_VAL Summary of this function goes here
%   Detailed explanation goes here

%% inner edge - default
% ind = (ftype == ndg_lib.bc_type.Inner);
% cP(ind) = cP(ind);

%% slip wall and non-slip wall
ind = (ftype == ndg_lib.bc_type.SlipWall);
cP(ind) = cM(ind);
ind = (ftype == ndg_lib.bc_type.SlipWall);
cP(ind) = cM(ind);

%% zero gradient
ind = (ftype == ndg_lib.bc_type.ZeroGrad);
cP(ind) = cM(ind);

%% clamped
ind = (ftype == ndg_lib.bc_type.Clamped);
cP(ind) = c_extM(ind);

end
