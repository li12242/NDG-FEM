function DrawSnapshot

% parameter
filename = 'Convection2D_1_40.nc';
sstime   = [2.4]; % snapshot time
eletype  = 'tri';
deg      = 1;

% dimensions and variables
time     = ncread(filename, 'time');

% read and draw
for t = 1:numel(sstime)
    [~, ist] = min(abs(time - sstime(t)));
    DrawPatch(filename, ist, eletype, deg);
    title(['$t=',num2str(sstime(t)),'$'], 'Interpreter', 'latex',...
        'FontSize', 18)
end% for

end% func

%% DrawPatch
% Draw snapshots at spicific time step
function DrawPatch(filename, ist, eleType, deg)
% dimensions and variables
varname  = 'var';
x        = ncread(filename, 'x');
y        = ncread(filename, 'y');

[np, ne] = size(x);

% get result
var      = ncread(filename, varname, [1,1,ist], [np, ne, 1]);

% draw result
figure;
vertex   = [x(:), y(:), var(:)];
EToV     = GetBcList(ne, deg, eleType);
patch('Vertices', vertex, 'Faces', EToV, 'FaceColor', [0.8, 0.9, 1])
axis equal
xlim([-1, 1]);
ylim([-1, 1]);
zlim([-0.2, 1.2]);
view(-29, 20)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
zlabel('$C$', 'Interpreter', 'latex')

end% func

%% GetBcList
% Return the boundary node list of each element
function EToV = GetBcList(ne, deg, eleType)
switch eleType
    case 'tri'
        tri    = StdRegions.Triangle(deg);
        bclist = tri.getFaceListToNodeList';
        EToV   = ones(ne, 1)*bclist;
        EToV   = EToV + [tri.nNode*(0:ne-1)]'*ones(size(bclist));
    case 'quad'
        quad   = StdRegions.Quad(deg);
        bclist = quad.getFaceListToNodeList;
        EToV   = ones(ne, 1)*bclist;
        EToV   = EToV + [tri.nNode*(0:ne-1)]'*ones(size(bclist));
end% switch
end% func