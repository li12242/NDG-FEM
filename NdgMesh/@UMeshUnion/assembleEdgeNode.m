function [eidM, eidP, eidtype] = assembleEdgeNode(obj)
TNfp = obj.cell.TNfp;
idm = zeros(TNfp, obj.K);
idp = zeros(TNfp, obj.K);
type = NdgEdgeType.Inner * ones(TNfp, obj.K, 'int8');

Np = obj.cell.Np;
for k1 = 1:obj.K
    idm(:, k1) = (k1-1)*Np + obj.cell.Fmask(:);
    f1 = 1; nfp1 = 1;
    for n = 1:TNfp
        if nfp1 > obj.cell.Nfp(f1)
            f1 = f1 + 1;
            nfp1 = 1;
        end
        type(n, k1) = obj.EToB(f1, k1);
        
        k2 = obj.EToE(f1, k1);
        f2 = obj.EToF(f1, k1);
        
        xpM = obj.x(idm(n, k1));
        ypM = obj.y(idm(n, k1));
        zpM = obj.z(idm(n, k1));
        ind2 = (k2-1)*Np + obj.cell.Fmask(:, f2);
        xP = obj.x( ind2 );
        yP = obj.y( ind2 );
        zP = obj.z( ind2 );
        d12 = (xpM - xP).^2 + (ypM - yP).^2 + (zpM - zP).^2;
        m = (d12 < 3e-16);
        try
        idp(n, k1) = (k2-1)*Np + obj.cell.Fmask(m, f2);
        catch 
            keyboard
        end
        nfp1 = nfp1 + 1;
    end
end
eidM = idm;
eidP = idp;
eidtype = type;

end