function [ fval ] = orthogonal_func(obj, N, td, r, s, t)
%ORTHOGONAL_FUNC Get the values of the orthgonal basis
%   Get the i-th orthgonal function value at the coordinate (r,s,t)

% project the coordinate (r,s) in triangle to (a,b) in square
[ a,b ] = rstoab( r,s ); 

% 基函数序号转换为矩形内正交函数编号（i，j）
[ i, j ] = trans_ind( N, td );
[ fval ] = simplex2DP( a,b,i,j );

end

function [ P ] = simplex2DP( a,b,i,j )
% Evaluate 2D orthonormal polynomial on simplex at (a,b) of order (i,j).
% 
h1 = JacobiP(a,0,0,i); 
h2 = JacobiP(b,2*i+1,0,j);
P = sqrt(2.0)*h1.*h2.*(1-b).^i;
end



