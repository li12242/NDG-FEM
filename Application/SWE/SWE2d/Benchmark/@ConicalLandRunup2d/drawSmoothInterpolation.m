sx = linspace( 0, 25, 200 );
sy = linspace( 0, 30, 200 );

[sx, sy] = meshgrid(sx, sy);

a1 = gca;
Faces = a1.Children(1).Faces;
Vertices = a1.Children(1).Vertices;

x = a1.Children(1).XData;
y = a1.Children(1).YData;
z = a1.Children(1).ZData;

F = TriScatteredInterp(x(:), y(:), z(:), 'nearest');
%sz = reshape( F(sx, sy), 200, 200 );

%Fs = TriScatteredInterp(sx, sy, sz);
z = F(Vertices(:,1), Vertices(:,2));
Vertices(:, 3) = z(:);

set(a1.Children(1), 'Vertices', Vertices,...
    'FaceVertexCData', z(:));

