function [ EToM, EToE, EToF, EToB ] = mxAssembleMeshConnection( ...
    MeshId1, MeshId2, K1, K2, ...
    Nface1, Nface2, FToV1, FToV2, EToV1, EToV2, ...
    EToM, EToE, EToF, EToB )
%MXASSEMBLEMESHCONNECTION Summary of this function goes here
%   Detailed explanation goes here

for k = 1 : K1 
    for f = 1 : Nface1
        EToM(f, k) = MeshId1;
        if EToE(f, k) == k
            v11 = max( EToV1(FToV1(1:2, f), k) );
            v12 = min( EToV1(FToV1(1:2, f), k) );

            for k2 = 1 : K2
                for f2 = 1 : Nface2
                    vid = EToV2(FToV2(1:2, f2), k2);
                    v21 = max( vid );
                    v22 = min( vid );

                    if (v21 == v11) && (v12 == v22)
                        EToM(f, k) = MeshId2;
                        EToE(f, k) = k2;
                        EToF(f, k) = f2;
                        EToB(f, k) = 1; % Gauss edge
                    end
                end
            end
        end
    end
end

end

