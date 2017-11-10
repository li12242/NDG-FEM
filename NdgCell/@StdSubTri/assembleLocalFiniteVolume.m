function [ NFV, EToV ] = assembleLocalFiniteVolume( obj )

NFV = obj.N^2;
EToV = zeros(3, NFV);

sk = 1;
for row = 1:obj.N
    v1 = obj.Fmask(row, 3);
    v2 = obj.Fmask(row+1, 3);
    for kb = 1:row
        EToV(:, sk) = [v1, v2, v2+1]';
        v1 = v1+1;
        v2 = v2+1;
        sk = sk+1;
    end
    
    v1 = obj.Fmask(row, 3);
    v2 = obj.Fmask(row+1, 3);
    for kt = 1:row-1
        EToV(:, sk) = [v1, v2+1, v1+1]';
        v1 = v1+1;
        v2 = v2+1;
        sk = sk+1;
    end
end
end

