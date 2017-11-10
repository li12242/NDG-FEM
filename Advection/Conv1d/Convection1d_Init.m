function phys = Convection1d_Init(phys)

switch phys.casename
    case 'Advection'
        var = Advection_Init(phys);
    case 'AdvectionDiffusion'
        var = AdvectionDiffusion_Init(phys);
end
phys.var = var;
end% func

function var = AdvectionDiffusion_Init(phys)
x1 = phys.xlim(1); 
x2 = phys.xlim(2); 
x0 = x1 + 0.75*(x2 - x1);

mesh = phys.mesh;
Dx = phys.Dx;
t = -(mesh.x-x0).^2/Dx;
var = exp(t);
end

function var = Advection_Init(phys)
x1 = phys.xlim(1); 
x2 = phys.xlim(2); 
x0 = (x1+x2)./2;
length = (x2 - x1);

mesh = phys.mesh;
var = zeros(size(mesh.x));
% left part
flag = (mesh.x <= x0); x11 = (x1 + x0)./2;
var(flag) = exp( - (mesh.x(flag) - x11).^2./0.005/length^2 );

% right part
% x21 = (x2 + 3*x0)/4; x22 = (3*x2 + x0)/4;
% flag = (mesh.x >= x21) & (mesh.x <= x22);
% var(flag) = 1;
end