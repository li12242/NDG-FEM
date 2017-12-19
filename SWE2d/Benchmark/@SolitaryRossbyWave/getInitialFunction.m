function [ fext ] = getInitialFunction( obj )
fext = cell( obj.Nmesh, 1 );

for m = 1:obj.Nmesh%网格循环
    fext{m} = zeros(obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield);%大小Np*K*Nfield
    a = obj.meshUnion(m).x;
    b = obj.meshUnion(m).y;
    [ h, u, v ] = getInitialRossbyWave(obj, obj.meshUnion(m), a, b);
    hu = h.*u;
    hv = h.*v;
    fext{m}(:,:,1) = h;
    fext{m}(:,:,2) = hu;
    fext{m}(:,:,3) = hv;
end

end%func

function [ h, u, v ] = getInitialRossbyWave(obj, mesh, x, y)
h     = ones(mesh.cell.Np, mesh.K);
u     = zeros(mesh.cell.Np, mesh.K);
v     = zeros(mesh.cell.Np, mesh.K);
s     = obj.a;

r     = 0.771*(s^2)*(sech(s*x)).*(sech(s*x));
dr    = 0.771*(-2)*(s^3)*(sech(s*x).*sech(s*x)).*tanh(s*x);

h = 1 + r.*((3+6*(y.*y))/4).*(exp(-(y.*y)/2));
u = r.*((-9+6*(y.*y))/4).*(exp(-(y.*y)/2));
v = dr.*(2*y).*(exp(-(y.*y)/2));


% r1     = 0.771*(s^2)*(sech(s*((3+6*(y.*y))/4))).*(sech(s*((3+6*(y.*y))/4)));
% r2     = 0.771*(s^2)*(sech(s*((-9+6*(y.*y))/4))).*(sech(s*((-9+6*(y.*y))/4)));
% dr    = 0.771*(-2)*(s^3)*(sech(s*2*y)).*(sech(s*2*y)).*tanh(s*2*y);
% 
% % h = 1 + r1.*(exp(-(y.*y)/2));
% u = r2.*(exp(-(y.*y)/2));
% v = dr.*(exp(-(y.*y)/2));


end%func

