classdef test_ugrid_mesh2d < matlab.unittest.TestCase
    methods (Test)
        function testTmpMaxBelowThreshold(testCase)
            % Setup
            setup();

            % Define the mesh file path
            meshfile = 'Application/SWE/SWE2d/Benchmark/@ConicalLandRunup2d/mesh/conicalLand.msh';
            mesh = read_gmsh2d(meshfile);

            % Create UGRID sets
            vert = ugrid_set(mesh.Nv, "vert");
            cell = ugrid_set(mesh.Ntri, "cell");

            % Define UGRID data
            vx = ugrid_data(vert, 1, ugrid_dataTypeEnum.Float, mesh.vx, "vx");
            xc = ugrid_data(cell, 1, ugrid_dataTypeEnum.Float, zeros(mesh.Ntri, 1), "xc");

            % Create mapping
            imap = sparse(repmat(1:mesh.Ntri, 3, 1), mesh.EToVTri, 1/3, ...
                mesh.Ntri, mesh.Nv);
            map_v2c = ugrid_map(vx, xc, 3*mesh.Ntri, imap, "map_v2c");
            xc.value = map_v2c.imap * vx.value;

            % Calculate tmp
            tmp = max(abs(xc.value - sum(mesh.vx(mesh.EToVTri'), 2)/3));

            % Verify tmp is below the threshold
            threshold = 1e-6; % Define the threshold value
            testCase.verifyLessThan(tmp, threshold, "tmp exceeds the threshold value");
        end
    end
end

% --------------------------------------------------------------
function setup()
addpath(genpath('src/ugrid/src/interface/matlab'));
end % function setup()

% --------------------------------------------------------------
function mesh = read_gmsh2d(meshfile)
% READ_GMSH2D Read a 2D mesh from a Gmsh file

% read mesh data with ugrid solver
fid1 = fopen(meshfile, 'r');
if( fid1 < 0 )
    msgID = [mfilename, ':inputFileNameError'];
    msgtext = ['The input file name: ', meshfile, ' is incorrect'];
    ME = MException(msgID, msgtext);
    throw(ME);
end

% Skip header lines
for i = 1:4, fgetl(fid1); end

% Read vertex data
Nv = fscanf(fid1, '%d', 1);
data = fscanf(fid1, '%*d %f %f %*f\n', [2, Nv]);
vx = data(1, :)';
vy = data(2, :)';

% Skip to element data
for i = 1:3, fgetl(fid1); end
Ne = fscanf(fid1, '%d\n', 1);

% Initialize counters
Nedge = 0; Ntri = 0; Nquad = 0;

% Count element types
for i = 1:Ne
    elemType = sscanf(fgetl(fid1), '%*d %d', 1);
    switch elemType
        case 1, Nedge = Nedge + 1;
        case 2, Ntri = Ntri + 1;
        case 3, Nquad = Nquad + 1;
    end
end

% % Output values for verification
% fprintf('Nv: %d\n', Nv);
% fprintf('Nedge: %d\n', Nedge);
% fprintf('Ntri: %d\n', Ntri);
% fprintf('Nquad: %d\n', Nquad);

fseek(fid1,0,-1);
for i=1:(Nv+8), fgetl(fid1); end

% data = fscanf(fid1,'%d %d %d %d %d %d %d\n',[7,Nedge]);
% BCToV = [data(6, :); data(7, :); data(4, :)];
data = fscanf(fid1, '%*d %*d %*d %d %*d %d %d\n', [3, Nedge]);
BCToV = [data(2:3, :); data(1, :)];

% fprintf("data size %d %d\n", size(data, 1), size(data, 2));
% fprintf("BCToV size %d %d\n", size(BCToV, 1), size(BCToV, 2));
% fprintf("BCToV sample %d %d %d\n", BCToV(:,end));

mesh = struct('Nv', Nv, 'vx', vx, 'vy', vy, 'Nedge', Nedge, 'BCToV', BCToV);

% Read triangle connectivity
if Ntri > 0
    data = fscanf(fid1, '%*d %*d %*d %d %*d %d %d %d\n', [4, Ntri]);
    mesh.Ntri = Ntri;
    mesh.EToVTri = data(2:4, :);
    mesh.regid_tri = data(1, :);
end

% Read quadrilateral connectivity
if Nquad > 0
    data = fscanf(fid1, '%*d %*d %*d %d %*d %d %d %d %d\n', [5, Nquad]);
    mesh.Nquad = Nquad;
    mesh.EToVQuad = data(2:5, :);
    mesh.regid_quad = data(1, :);
end

end % function read_gmsh2d