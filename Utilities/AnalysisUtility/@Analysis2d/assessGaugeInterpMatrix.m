function [ Vg ] = assessGaugeInterpMatrix( obj )
    Vg = cell(obj.Ng, 1);

    for n = 1 : obj.Ng
        cellId = obj.gaugeCell(n);
        m = obj.gaugeMesh(n);
        mesh = obj.phys.meshUnion(m);

        if (mesh.cell.type == enumStdCell.Tri)
            [rd, sd] = accessTriLocalCoor( mesh, obj.xg(n), obj.yg(n), cellId );
        elseif (mesh.cell.type == enumStdCell.Quad)
            [rd, sd] = accessQuadLocalCoor( mesh, obj.xg(n), obj.yg(n), cellId );
        end

        Vg{n} = mesh.cell.nodal_func( rd, sd, 0 );
    end

end

function [rd, sd] = accessQuadLocalCoor( mesh, xd, yd, kd )

fmask1 = mesh.cell.Fmask(1, :)';
fmask2 = mesh.cell.Fmask(end, :)';

xv1 = mesh.x( fmask1, kd );
yv1 = mesh.y( fmask1, kd );
xv2 = mesh.x( fmask2, kd );
yv2 = mesh.y( fmask2, kd );
% obtain the local coordinate by the area coordinate
b = sqrt( (xv1 - xv2).^2 + (yv1 - yv2).^2 );
a1 = sqrt( (xd - xv1).^2 + (yd - yv1).^2 );
a2 = sqrt( (xd - xv2).^2 + (yd - yv2).^2 );
d = ( b + a1 + a2 )./2;
L = sqrt( d.*(d-b).*(d-a1).*(d-a2) )./mesh.LAV( kd );
rd = (L(4) - L(2))*2;
sd = (L(1) - L(3))*2;

end

function [rd, sd] = accessTriLocalCoor( mesh, xd, yd, kd )

fmask1 = mesh.cell.Fmask(1, :)';
fmask2 = mesh.cell.Fmask(end, :)';

xv1 = mesh.x( fmask1, kd );
yv1 = mesh.y( fmask1, kd );
xv2 = mesh.x( fmask2, kd );
yv2 = mesh.y( fmask2, kd );
% obtain the local coordinate by the area coordinate
b = sqrt( (xv1 - xv2).^2 + (yv1 - yv2).^2 );
a1 = sqrt( (xd - xv1).^2 + (yd - yv1).^2 );
a2 = sqrt( (xd - xv2).^2 + (yd - yv2).^2 );
d = ( b + a1 + a2 )./2;
L = sqrt( d.*(d-b).*(d-a1).*(d-a2) )./mesh.LAV( kd );
rd = -L(2) + L(3) - L(1);
sd = -L(2) - L(3) + L(1);
end