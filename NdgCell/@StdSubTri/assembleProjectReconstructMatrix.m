function [ P, R ] = assembleProjectReconstructMatrix( obj )

Nface = obj.Nface;
r = obj.r; 
s = obj.s;
P = zeros( obj.Np, obj.Np );
Pt = zeros( obj.Nq, obj.Np);

for k = 1:obj.NFV
    vertId = obj.EToV(:, k);
    rc = mean( r(vertId) );
    sc = mean( s(vertId) );
    
    fvArea = TriArea( r(vertId), s(vertId) );
    for f1 = 1:Nface
        vind1 = f1;
        vind2 = mod(vind1, obj.Nv)+1;
        v1 = vertId( vind1 );
        v2 = vertId( vind2 );
        
        r1 = r(v1); r2 = r(v2); rf = 0.5*(r1+r2);
        s1 = s(v1); s2 = s(v2); sf = 0.5*(s1+s2);
        
        rq1 = obj.project_vert2quad([rf, rc, r1]');
        sq1 = obj.project_vert2quad([sf, sc, s1]');
        
        for i = 1:obj.Np
            Pt(:,i) = obj.orthogonal_func(obj.N, i, rq1, sq1);
        end
        P(v1, :) = P(v1, :) + (fvArea/6/2 .*obj.wq')*Pt;
        
        rq2 = obj.project_vert2quad([rf, r2, rc]');
        sq2 = obj.project_vert2quad([sf, s2, sc]');
        for i = 1:obj.Np
            Pt(:,i) = obj.orthogonal_func(obj.N, i, rq2, sq2);
        end
        P(v2, :) = P(v2, :) + (fvArea/6/2 .*obj.wq')*Pt;
    end
end
P = P/obj.V;
P = bsxfun(@times, P, 1./obj.LAVToCV);
R = inv(P);
end

