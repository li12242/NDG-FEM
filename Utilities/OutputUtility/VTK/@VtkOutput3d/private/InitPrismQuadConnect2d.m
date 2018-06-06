function [ Ncell, EToV ] = InitQuadConnect2d(N, Nz)
    %InitQuadConnect2d - Description
    %
    % Syntax: [ Ncell, EToV ] = InitQuadConnect2d(input)
    %
    % Long description
    Np = ( N + 1 ) .^ 2;
    Ncell = N ^ 2 * Nz;
    EToV = zeros(8, Ncell);
    
    sk = 1; 
    for lay = 1 : Nz
        s2 = Np - N + Np * (lay - 1);
        s1 = s2 - ( N + 1 ) + Np * (lay - 1);
        for row = 1 : N
            v1 = s1; v2 = s2;
            for kb = 1:N
                EToV(:, sk) = [ v1, v1 + 1, v2+1, v2, ....
                    v1 + Np, v1 + Np + 1, v2 + Np + 1, v2 + Np ]';
                v1 = v1+1; v2 = v2+1; sk = sk+1;
            end
        
            s1 = s1 - (N + 1);
            s2 = s2 - (N + 1);
        end
    end
end