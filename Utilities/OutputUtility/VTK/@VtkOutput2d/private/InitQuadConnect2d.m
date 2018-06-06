function [ Ncell, EToV ] = InitQuadConnect2d(N)
%InitQuadConnect2d - Description
%
% Syntax: [ Ncell, EToV ] = InitQuadConnect2d(input)
%
% Long description
Np = ( N + 1 ) .^ 2;
Ncell = N ^ 2;
EToV = zeros(4, Ncell);

sk = 1; 
s2 = Np - N; 
s1 = s2 - ( N + 1 );
for row = 1 : N
    v1 = s1; v2 = s2;
    for kb = 1:N
        EToV(:, sk) = [ v1, v1 + 1, v2+1, v2 ]';
        v1 = v1+1; v2 = v2+1; sk = sk+1;
    end

    s1 = s1 - (N + 1);
    s2 = s2 - (N + 1);
end
end