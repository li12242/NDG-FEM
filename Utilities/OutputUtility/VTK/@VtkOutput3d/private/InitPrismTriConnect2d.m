function [ Ncell, EToV ] = InitPrismTriConnect2d( N, Nz )
    
    Np = ( N + 1 ) * ( N + 2 ) / 2;
    Ncell = N ^ 2 * Nz;
    EToV = zeros(6, Ncell);

    sk = 1; 
    for lay = 1 : Nz
        s2 = Np + Np * (lay - 1);
        s1 = Np - 2 + Np * (lay - 1);
        for row = 1 : N
            v1 = s2; v2 = s1;
            for kb = 1:row
                EToV(:, sk) = [v1, v2, v2 + 1, ...
                    v1 + Np, v2 + Np, v2 + Np + 1]';
                v1 = v1+1; v2 = v2+1; sk = sk+1;
            end
            
            v1 = s2; v2 = s1;
            for kt = 1:row-1
                EToV(:, sk) = [v1, v2 + 1, v1 + 1, ...
                    v1 + Np, v2 + Np + 1, v1 + Np + 1]';
                v1 = v1+1; v2 = v2+1; sk = sk+1;
            end

            s1 = s1 - (row + 2);
            s2 = s2 - (row + 1);
        end
    end
end