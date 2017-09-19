function [ dr, ds ] = tri_deri_orthgonal_func( i, j, r, s )
%TRI_DERI_ORTHGONAL_FUNC Summary of this function goes here
%   Detailed explanation goes here

% project the coordinate from the triangle to the square
[ a,b ] = rstoab(r,s);
% 计算对（r,s）坐标的导数
[dr, ds] = deri_simplex2DP(a,b,i,j);
end

function [ dfdr, dfds ] = deri_simplex2DP( a,b,id,jd )
% Return the derivatives of the modal basis (id,jd) on the 2D simplex at (a,b).
fa = Polylib.JacobiP(a, 0, 0, id);     
dfa = Polylib.GradJacobiP(a, 0, 0, id);
gb = Polylib.JacobiP(b, 2*id+1,0, jd); 
dgb = Polylib.GradJacobiP(b, 2*id+1,0, jd);
% r-derivative
% d/dr = da/dr d/da + db/dr d/db = (2/(1-s)) d/da = (2/(1-b)) d/da
dfdr = dfa.*gb;
if(id>0)
    dfdr = dfdr.*((0.5*(1-b)).^(id-1));
end
% s-derivative
% d/ds = ((1+a)/2)/((1-b)/2) d/da + d/db
dfds = dfa.*(gb.*(0.5*(1+a)));
if(id>0)
    dfds = dfds.*((0.5*(1-b)).^(id-1));
end
tmp = dgb.*((0.5*(1-b)).^id);
if(id>0)
    tmp = tmp-0.5*id*gb.*((0.5*(1-b)).^(id-1));
end
dfds = dfds+fa.*tmp;
% Normalize
dfdr = 2^(id+0.5)*dfdr; dfds = 2^(id+0.5)*dfds;
end

