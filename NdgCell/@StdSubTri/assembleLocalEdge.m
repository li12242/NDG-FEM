function [ Nedge, nx, ny, v1, v2 ] = assembleLocalEdge( obj )

Nface = obj.Nface;
Nedge = obj.NFV*Nface;
v1 = zeros(Nface, obj.NFV);
v2 = zeros(Nface, obj.NFV);
nx = zeros(Nface, obj.NFV);
ny = zeros(Nface, obj.NFV);
for k = 1:obj.NFV % loop over FVs
    vertId = obj.EToV(:, k);
    xc = mean( obj.r( vertId ) );
    yc = mean( obj.s( vertId ) );
    
    for f1 = 1:Nface
        f2 = mod(f1,Nface)+1;
        fv1 = vertId( f1 );
        fv2 = vertId( f2 );
        
        v1(f1, k) = fv1;
        v2(f1, k) = fv2;
        
        x1 = obj.r( fv1 );
        x2 = obj.r( fv2 );
        y1 = obj.s( fv1 );
        y2 = obj.s( fv2 );
        
        xf = 0.5*( x1 + x2 );
        yf = 0.5*( y1 + y2 );
        nxt = yc - yf;
        nyt = xf - xc;
        
        % check the direction is from v1 to v2
        dx = x2 - x1;
        dy = y2 - y1;
        ds = sqrt( nxt.^2 + nyt.^2 );
        if ( dx*nxt + dy*nyt ) < 0
            nx(f1, k) = -nxt/ds;
            ny(f1, k) = -nyt/ds;
        else
            nx(f1, k) = nxt/ds;
            ny(f1, k) = nyt/ds;
        end
    end
end
v1 = v1(:);
v2 = v2(:);

nx = nx(:);
ny = ny(:);
end

