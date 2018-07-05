function [ f ] = EvaluateVerticalOrthogonalFunc( obj, td, t )
%EVALUATEVERTICALORTHOGONALFUNC Summary of this function goes here
%   Detailed explanation goes here

td2 = ceil( td / obj.Nph );
[ f ] = line_orthogonal_func( td2, t );

end

function [ f ] = line_orthogonal_func( ind, r )
f = JacobiP(r, 0, 0, ind-1);
end% func
