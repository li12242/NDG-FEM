function [ Ncell, EToV ] = InitTriConnect2d(N)
%InitTriConnect2d - Description
%
% Syntax: [ Ncell, EToV ] = InitTriConnect2d(mesh)
%
% Long description
    
    Np = ( N + 1 ) * ( N + 2 ) / 2;
    Ncell = N ^ 2;
    EToV = zeros(3, Ncell);

    sk = 1; 
    s2 = Np; 
    s1 = Np - 2;
    for row = 1 : N
        v1 = s2; v2 = s1;
        for kb = 1:row
            EToV(:, sk) = [v1, v2, v2+1]';
            v1 = v1+1; v2 = v2+1; sk = sk+1;
        end
        
        v1 = s2; v2 = s1;
        for kt = 1:row-1
            EToV(:, sk) = [v1, v2+1, v1+1]';
            v1 = v1+1; v2 = v2+1; sk = sk+1;
        end

        s1 = s1 - (row + 2);
        s2 = s2 - (row + 1);
    end

end