function obj = ConnectMesh( obj )

Nface = obj.cell.Nface;
Ke = obj.K;
Nfv = max( obj.cell.Nfv );
faceVertIndex = zeros( Nfv + 2, Nface * Ke );

sk = 1;
for k = 1 : Ke
    for f = 1 : Nface
        faceVertIndex( 1:Nfv, sk ) = sort( obj.EToV( obj.cell.FToV(:, f) ,k) );
        faceVertIndex( Nfv + 1, sk ) = k;
        faceVertIndex( Nfv + 2, sk ) = f;
        sk = sk + 1;
    end
end

faceVertIndex = faceVertIndex';
faceVertIndex = sortrows( faceVertIndex );
faceVertIndex = faceVertIndex';

EToE = ones(Nface, 1) * (1 : Ke);
EToF = (1 : Nface)' * ones(1, Ke);

for i = 1 : ( Ke * Nface - 1 )
    if all( faceVertIndex(1 : Nfv, i + 1) == faceVertIndex(1 : Nfv, i) )
        k1 = faceVertIndex( Nfv + 1, i );
        f1 = faceVertIndex( Nfv + 2, i );
        k2 = faceVertIndex( Nfv + 1, i + 1 );
        f2 = faceVertIndex( Nfv + 2, i + 1 );
        EToE(f1, k1) = k2; 
        EToE(f2, k2) = k1;
        EToF(f1, k1) = f2; 
        EToF(f2, k2) = f1;
    end
end

obj.EToE = EToE;
obj.EToF = EToF;
end
