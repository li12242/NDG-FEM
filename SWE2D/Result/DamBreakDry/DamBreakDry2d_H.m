%% DamBreakQx
% Exact solution for water depth
function h = DamBreakDry2d_H(x, y, t)
h0    = 10;
xc    = 500;
theta = (x - xc)/t;
g     = 9.81; 
sgh   = sqrt(g*h0);

h     = zeros(size(x));

% wet part
ind    = theta < -sgh;
h(ind) = h0;

ind    = theta > 2*sgh;
h(ind) = 0;

ind    = (theta >= -sgh) & (theta <= 2*sgh);
h(ind) = 1/9/g*(theta(ind) - 2*sgh).^2;

end% func