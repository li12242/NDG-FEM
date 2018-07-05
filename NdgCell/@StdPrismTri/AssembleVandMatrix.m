function AssembleVandMatrix( obj )
V = zeros(obj.Np, obj.Np);

for n = 1:obj.Np
    fh = obj.EvaluateHorizontalOrthogonalFunc( obj.N, n, obj.r, obj.s );
    fv = obj.EvaluateVerticalOrthogonalFunc( n, obj.t );
    V(:, n) = fh .* fv;
end% for

Vh = zeros(obj.Nph, obj.Nph);
for n = 1:obj.Nph
    fh = obj.EvaluateHorizontalOrthogonalFunc( obj.N, n, obj.r1, obj.s1 );
    Vh(:, n) = fh;
end% for

% vertical integral vandermonde matrix
Vint = EvaluateVerticalIntegralOrthogonalFunc( obj );

obj.V = V;
obj.Vh = Vh;
obj.Vint = Vint;

end% func

function Vz = EvaluateVerticalIntegralOrthogonalFunc( obj )

Vz = zeros(obj.Np, obj.Np);
for n = 1:obj.Np
    fh = obj.EvaluateHorizontalOrthogonalFunc( obj.N, n, obj.r, obj.s );
    ind = ceil( n / obj.Nph ); % vertical orthogonal polynomial index
    if ind == 1
        fv = LegendreNorm1d( ind + 1, obj.t );
        Vz(:, n) = ( fv - fv(1) ) ./ sqrt( 2 * ind + 1 );
    else
        fvm = LegendreNorm1d( ind - 1, obj.t );
        fvp = LegendreNorm1d( ind + 1, obj.t );
        Vz(:, n) = ( fvp ./ sqrt( 2*ind+1 ) - fvm ./ sqrt( 2*ind-3 ) ) ./ sqrt( 2*ind-1 );
    end
    
    Vz(:, n) = Vz(:, n) .* fh;
end% for

end

function [ f ] = LegendreNorm1d( ind, r )
f = JacobiP(r, 0, 0, ind-1);
end% func