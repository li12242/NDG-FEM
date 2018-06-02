adv = ConstAdvUniformMesh2d( 1, 4, NdgCellType.Tri );
mesh = adv.meshUnion;
std = StdPrismTri( 1, 8 );
zs = zeros(mesh.Nv, 1);
zb = zeros(mesh.Nv, 1) - 1;
mesh3 = NdgExtendMesh3d( std, mesh, zs, zb, 2 );
edge3 = NdgSideEdge3d( mesh3, 1 );
bEdge3 = NdgBottomEdge3d( mesh3, 1 );

xc = - 0.5; yc = 0; zc = -0.5;
r2 = sqrt( (mesh3.x - xc).^2 + (mesh3.y - yc).^2 ...
    + (mesh3.z - zc).^2 ) ./ 0.25;
ind = ( r2 <= 1.0);
f = cell(1, 1);
f{1} = zeros( mesh3.cell.Np, mesh3.K );
f{1}(ind) = ( 1+cos( r2(ind)*pi ) )./2;

% scatter3( mesh3.x(:), mesh3.y(:), mesh3.z(:), f(:)*20 + 2, f(:)*20 + 2 )
Nlayer = 6;
ind = [1,2,3] + 3 * (Nlayer - 1);
mesh.draw( f{1}(ind, :) )

rhs = zeros( mesh3.cell.Np, mesh3.K );
u = 0.5;
dt = 0.5 * min( mesh.LAV / u );
Nt = 1 / (dt * u );


for i = 1 : Nt
    [ fm, fp ] = edge3.matEvaluateSurfValue( f );
    fluxM = u .* fm .* edge3.nx;
    fluxP = u .* fp .* edge3.nx;
    % numerical flux
    uNorm = u .* edge3.nx;
    sign_um = sign( uNorm );
    fluxS = ( fm .* ( sign_um + 1 ) .* 0.5 + fp .* ( 1 - sign_um  ) .* 0.5 ) .* uNorm;
    rhs = edge3.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS );
    
    E = u .* f{1};
    rhs = rhs - mesh3.rx.*( mesh3.cell.Dr * E ) ...
        - mesh3.sx.*( mesh3.cell.Ds * E ) ...
        - mesh3.tx.*( mesh3.cell.Dt * E );
    
    f{1} = f{1} + dt * rhs;
end
