function [ cellId, Vg ] = accessGaugePointLocation( obj, xg, yg, zg )

Ng = numel( xg );
xb = obj.x([1, end], :);
rd = zeros( Ng, 1 );
cellId = zeros( Ng, 1 );
for n = 1:Ng
    dx = xb - xg(n);
    id = find( dx(1,:).*dx(2,:) <= 0, 1 );
    cellId(n) = id;
    rd(n) = ( xg(n) - obj.x(1, id) )./( obj.x(2, id) - obj.x(1, id) )*obj.cell.LAV - 1;
end

Vg = zeros(Ng, obj.cell.Np);
for n = 1:obj.cell.Np
    Vg(:, n) = obj.cell.orthogonal_func(obj.cell.N, n, rd, 0, 0);
end
Vg = Vg/obj.cell.V;

end

