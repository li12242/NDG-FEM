function [ f ] = evaluateVerticalOrthogonalFunc( obj, td, t )
%EVALUATEVERTICALORTHOGONALFUNC Summary of this function goes here
%   Detailed explanation goes here

td2 = ceil( td / obj.Nph );
[ f ] = evaluateLineOrthogonalFunc( td2, t );

end

function [ f ] = evaluateLineOrthogonalFunc( ind, r )
f = JacobiP(r, 0, 0, ind-1);
end% func
