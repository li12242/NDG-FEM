function var = Convection1DInit(mesh, x1, x2)
x0 = (x1+x2)./2;
length = (x2 - x1);

var = zeros(size(mesh.x));

% left part
flag = (mesh.x <= x0); x11 = (x1 + x0)./2;
var(flag) = exp( - (mesh.x(flag) - x11).^2./0.005/length^2 );

% right part
x21 = (x2 + 3*x0)/4; x22 = (3*x2 + x0)/4;
flag = (mesh.x >= x21) & (mesh.x <= x22);
var(flag) = 1;

end% func