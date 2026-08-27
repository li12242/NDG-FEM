function [ Np, r, s, t ] = node_coor_func( obj, nOrder )
%> @brief Distribution of the interpolation nodes on the reference triangle.
%> The shared warp & blend placement is in ndgcell.triNodeCoor.
[ Np, r, s ] = ndgcell.triNodeCoor( nOrder );
t = zeros(size(r));
end
