function [dr, ds, dt] = derivative_orthogonal_func( obj, N1, N2, td, r, s, t)

    td1 = mod( td-1, obj.Nph ) + 1;
    td2 = ceil( td / obj.Nph );

    [ dr, ds ] = tri_derivative_orthogonal_func( N1, td1, r, s);
    [ ft ] = line_orthogonal_func( td2, t );
    % multiply the vertical polynomial
    dr = dr .* ft;
    ds = ds .* ft;

    [ dt ] = line_derivative_orthogonal_func( td2, t);
    [ frs ] = tri_orthogonal_func( N1, td1, r, s );
    % multiply the horizontal polynomial
    dt = dt .* frs;
end

function [ f ] = line_orthogonal_func( ind, r )
    f = JacobiP(r, 0, 0, ind-1);
end% func

function [ fval ] = tri_orthogonal_func( N, td, r, s )
%ORTHOGONAL_FUNC Get the values of the orthgonal basis
%   Get the i-th orthgonal function value at the coordinate (r,s,t)

    % project the coordinate (r,s) in triangle to (a,b) in square
    [ a,b ] = rstoab( r,s ); 

    % �������ת��Ϊ�����������ţ�i��j��
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

function [dr] = line_derivative_orthogonal_func( ind, r )
    dr = GradJacobiP(r, 0, 0, ind-1);
end

function [dr, ds] = tri_derivative_orthogonal_func( N, ind, r, s )
%GRADORTHOGONALFUN

% project to square.
[a,b] = rstoab(r,s);

% transform the index to two indexes.
[i, j] = trans_ind(N,ind);

% calculate the derivative function vales.
[dr, ds] = deri_simplex2DP(a,b,i,j);
end

function [ dmodedr, dmodeds ] = deri_simplex2DP( a,b,id,jd )
% Return the derivatives of the modal basis (id,jd) on the 2D simplex at (a,b).
fa = JacobiP(a, 0, 0, id);     
dfa = GradJacobiP(a, 0, 0, id);
gb = JacobiP(b, 2*id+1,0, jd);
dgb = GradJacobiP(b, 2*id+1,0, jd);
% r-derivative
% d/dr = da/dr d/da + db/dr d/db = (2/(1-s)) d/da = (2/(1-b)) d/da
dmodedr = dfa.*gb;
if(id>0)
    dmodedr = dmodedr.*((0.5*(1-b)).^(id-1));
end
% s-derivative
% d/ds = ((1+a)/2)/((1-b)/2) d/da + d/db
dmodeds = dfa.*(gb.*(0.5*(1+a)));
if(id>0)
    dmodeds = dmodeds.*((0.5*(1-b)).^(id-1));
end
tmp = dgb.*((0.5*(1-b)).^id);
if(id>0)
    tmp = tmp-0.5*id*gb.*((0.5*(1-b)).^(id-1));
end
dmodeds = dmodeds+fa.*tmp;
% Normalize
dmodedr = 2^(id+0.5)*dmodedr; dmodeds = 2^(id+0.5)*dmodeds;
end

