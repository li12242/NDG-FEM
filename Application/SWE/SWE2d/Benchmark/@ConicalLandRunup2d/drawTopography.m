function drawTopography( obj )

M = 80;
x = linspace(0, 25, M);
y = linspace(0, 30, M);
[x, y] = meshgrid(x, y);

% topography
xt = 12.5; yt = 15;
r = sqrt( (x - xt).^2 + (y - yt).^2 );
bot = min( 0.625, 0.9 - r/4 );
bot( bot <= 0 ) = 0.0;

% colormap
figure('color', 'w')
cmap = colormap('jet');
% colormap(map)
contour(x, y, bot, 10, 'LineWidth', 1.5); hold on;
colorbar('northoutside'); 
xlabel('$x \rm{(m)}$', 'Interpreter', 'latex','FontSize', 16);
ylabel('$y \rm{(m)}$', 'Interpreter', 'latex','FontSize', 16);


gx = [ 6.36, 8.9, 9.9, 12.5, 15.1 ];            
gy = [ 14.25, 15, 15, 12.42, 15 ];
plot(gx, gy, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
axis equal;
box on;
grid on;
end

