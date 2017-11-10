function area = TriArea(vx, vy)
% calculate the area of triangle from given vertex coordinate
a = sqrt( (vx(1,:) - vx(2,:)).^2 + (vy(1,:)-vy(2,:)).^2 );
b = sqrt( (vx(2,:) - vx(3,:)).^2 + (vy(2,:)-vy(3,:)).^2 );
c = sqrt( (vx(3,:) - vx(1,:)).^2 + (vy(3,:)-vy(1,:)).^2 );
p = (a+b+c)./2;
area = sqrt(p.*(p-a).*(p-b).*(p-c));
end