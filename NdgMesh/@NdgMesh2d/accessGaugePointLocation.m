function [ cellId, Vg ] = accessGaugePointLocation( obj, xg, yg, zg )

Ng = numel( xg );
cellId = zeros( Ng, 1 );
Nv = obj.cell.Nv;
K = obj.K;
for n = 1:Ng
    vx = obj.vx( obj.EToV );
    vy = obj.vy( obj.EToV );
    
    gvx = xg(n) - vx;
    gvy = yg(n) - vy;
    
    cross = zeros( Nv, K );
    dot = zeros( Nv, K );
    for f1 = 1:Nv
        f2 = mod(f1, Nv) + 1;
        dvx = vx(f2, :) - vx(f1, :);
        dvy = vy(f2, :) - vy(f1, :);
        dot(f1, :) = ( gvx(f1, :) .* gvx(f1, :) + gvy(f1,:) .* gvy(f1,:) ) ...
            ./ ( dvx .* dvx + dvy .* dvy );
        cross(f1, :) = gvx(f1, :) .* gvy(f2, :) - gvx(f2, :) .* gvy(f1, :);
    end

    ind = max( cross ) .* min( cross );
    [ ~, k2 ] = find( ind >= 0, 1 );
    if ~isempty(k2)
        cellId( n ) = k2;
    end
end

% calculate the map matrix for each gauge points
Vg = zeros(Ng, obj.cell.Np);
switch obj.cell.type
    case NdgCellType.Tri
        [rd, sd] = accessTriLocalCoor( obj, xg, yg, cellId );
    case NdgCellType.Quad
        [rd, sd] = accessQuadLocalCoor( obj, xg, yg, cellId );
end

for n = 1:obj.cell.Np
    Vg(:, n) = obj.cell.orthogonal_func(obj.cell.N, n, rd, sd, 0);
end
Vg = Vg/obj.cell.V;
end

function [rd, sd] = accessQuadLocalCoor( obj, xd, yd, kd )
Ng = numel( xd );
rd = ones(Ng, 1);
sd = ones(Ng, 1);
fmask1 = obj.cell.Fmask(1, :)';
fmask2 = obj.cell.Fmask(end, :)';
for n = 1:Ng
    if (kd(n) == 0) 
        continue; 
    end
    xv1 = obj.x( fmask1, kd(n) );
    yv1 = obj.y( fmask1, kd(n) );
    xv2 = obj.x( fmask2, kd(n) );
    yv2 = obj.y( fmask2, kd(n) );
    % obtain the local coordinate by the area coordinate
    b = sqrt( (xv1 - xv2).^2 + (yv1 - yv2).^2 );
    a1 = sqrt( (xd(n) - xv1).^2 + (yd(n) - yv1).^2 );
    a2 = sqrt( (xd(n) - xv2).^2 + (yd(n) - yv2).^2 );
    d = ( b + a1 + a2 )./2;
    L = sqrt( d.*(d-b).*(d-a1).*(d-a2) )./obj.LAV( kd(n) );
    rd(n) = (L(4) - L(2))*2;
    sd(n) = (L(1) - L(3))*2;
end
end

function [rd, sd] = accessTriLocalCoor( obj, xd, yd, kd )
Ng = numel( xd );
rd = ones(Ng, 1);
sd = ones(Ng, 1);

fmask1 = obj.cell.Fmask(1, :)';
fmask2 = obj.cell.Fmask(end, :)';
for n = 1:Ng
    if (kd(n) == 0) 
        continue; 
    end
    xv1 = obj.x( fmask1, kd(n) );
    yv1 = obj.y( fmask1, kd(n) );
    xv2 = obj.x( fmask2, kd(n) );
    yv2 = obj.y( fmask2, kd(n) );
    % obtain the local coordinate by the area coordinate
    b = sqrt( (xv1 - xv2).^2 + (yv1 - yv2).^2 );
    a1 = sqrt( (xd(n) - xv1).^2 + (yd(n) - yv1).^2 );
    a2 = sqrt( (xd(n) - xv2).^2 + (yd(n) - yv2).^2 );
    d = ( b + a1 + a2 )./2;
    L = sqrt( d.*(d-b).*(d-a1).*(d-a2) )./obj.LAV( kd(n) );
    rd(n) = -L(2) + L(3) - L(1);
    sd(n) = -L(2) - L(3) + L(1);
end
end
