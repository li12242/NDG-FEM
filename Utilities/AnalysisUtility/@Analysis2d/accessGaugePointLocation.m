function [ gaugeMesh, gaugeCell ] = accessGaugePointLocation( obj )

    gaugeMesh = zeros( obj.Ng, 1 );
    gaugeCell = zeros( obj.Ng, 1 );
    for m = 1 : obj.phys.Nmesh
        mesh = obj.phys.meshUnion(m);
        Nv = mesh.cell.Nv;
        K = mesh.K;

        vx = mesh.vx( mesh.EToV );
        vy = mesh.vy( mesh.EToV );
        for n = 1 : obj.Ng
            
            gvx = obj.xg(n) - vx;
            gvy = obj.yg(n) - vy;
            
            cross = zeros( Nv, K );
            dot = zeros( Nv, K );
            for f1 = 1:Nv
                f2 = mod(f1, Nv) + 1;
                dvx = vx(f2, :) - vx(f1, :);
                dvy = vy(f2, :) - vy(f1, :);
                dot(f1, :) = ( gvx(f1, :) .* gvx(f1, :) ...
                    + gvy(f1,:) .* gvy(f1,:) ) ...
                    ./ ( dvx .* dvx + dvy .* dvy );
                cross(f1, :) = gvx(f1, :) .* gvy(f2, :) ...
                    - gvx(f2, :) .* gvy(f1, :);
            end

            ind = max( cross ) .* min( cross );
            [ ~, k2 ] = find( ind >= 0, 1 );
            if ~isempty(k2)
                gaugeMesh( n ) = m;
                gaugeCell( n ) = k2;
            end
        end
    end

end