%======================================================================
%> @brief Calculate the quadrature points and their weights in the standard
%> quadrilateral element.
%>
%> More detailed description.
%>
%> @param arg1 First argument
%> @param arg2 Second argument
%>
%> @retval out1 return value for the first output variable
%> @retval out2 return value for the second output variable
%======================================================================
%> This class is part of the NDG-FEM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
%======================================================================
function [Nq, rq, sq, tq, wq] = quad_coor_func(obj, N)
%QUADRATURE_NODE_FUNC Summary of this function goes here
%   Detailed explanation goes here

np = N+1;
% [zq, w] = zwgl(np); % the 1D LGL quadrature points and their weights
[ zq, w ] = zwglj(np);
% loop along the r-axis first
rq = zq*ones(1, np);
sq = ones(np, 1)*zq';
rq = rq(:); 
sq = sq(:); 
wq = bsxfun(@times, w, w'); 
wq = wq(:);
tq = zeros(size(rq));
Nq = numel(rq);
end

