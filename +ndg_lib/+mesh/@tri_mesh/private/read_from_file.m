function [Nv, vx, vy, K, EToV, EToR, EToBS] = read_from_file(casename)
%READ_FROM_FILE Summary of this function goes here
%   Detailed explanation goes here

isdebug = 0;
% read node file
[Nv, vx, vy] = read_node_file(casename);
vz = zeros(size(vx));
if isdebug
    fprintf('\n%s\n\n', 'In function read_from_file')
    fprintf('%s\n-- %s\n', 'finish read node file:', [casename, '.node']);
    fprintf('%s %d\n\n','Nv = ', Nv);
end
% read mesh file
[K, EToV, EToR] = read_mesh_file(casename);
EToV = resort_vert(EToV, vx, vy); % get the vertex anti-clockwise
if isdebug
    fprintf('%s\n-- %s\n', 'finish read element file:', [casename, '.ele']);
    fprintf('%s %d\n\n', 'K = ', K);
end
% read boundary condition file
EToBS = read_edge_file(casename, EToV, Nv);
if isdebug
    fprintf('%s\n-- %s\n\n', 'finish read boundary file:', [casename, '.edge']);
end
end% func

function [Nv, vx, vy] = read_node_file(casename)
filename = [casename, '.node'];
fp = fopen(filename);
% read vertex number Nv
Nv = fscanf(fp, '%d', 1);
fgetl(fp); % pass the rest of first line
% read vertex coordinate
data = fscanf(fp, '%d %g %g', [3, Nv]);
vx = data(2,:)';
vy = data(3,:)';

fclose(fp);
end% func

function [K, EToV, EToR] = read_mesh_file(casename)
filename = [casename, '.ele'];
fp = fopen(filename);
% read element number K
K = fscanf(fp, '%d', 1);
Nv = fscanf(fp, '%d', 1); % number of vertex in each cell
fgetl(fp); % pass the rest of first line
% read vertex list of each element
if Nv == 3
    fmatStr = '%d %d %d %d %d';
elseif Nv == 4
    fmatStr = '%d %d %d %d %d %d';
end
data = fscanf(fp, fmatStr, [Nv+2, K]);
EToV = data( 2:(1+Nv), :);
EToR = data( Nv+2, :);
fclose(fp);
end% func

function EToBS = read_edge_file(casename, EToV, Nv)
filename = [casename, '.edge'];
fp = fopen(filename);
% read total number of 
Nf = fscanf(fp, '%d', 1);
fgetl(fp); % pass the rest of first line

% read face to vertex list
data = fscanf(fp, '%d %d %d %d', [4, Nf]);
ind = min( data([2,3], :) )*Nv + max( data([2,3], :) );
ftype = data(4, :);

% assigned to EToBS
[Nface, K] = size(EToV);
EToBS = int8(ones(Nface, K)).*ndg_lib.bc_type.Inner; % initialize
for k = 1:K
    for f = 1:Nface
        v1 = EToV(f, k);
        v2 = EToV(mod(f, Nface)+1, k);
        t = min(v1, v2)*Nv + max(v1, v2);
        tnd = find( abs(ind - t)<1e-10 );
        if isempty(tnd)
            continue
        else
            EToBS(f, k) = ftype(tnd);
        end
    end
end

fclose(fp);
end% func