function [ obcfile ] = make_obc_file( obj, casename )
%MAKE_OBC_FILE Summary of this function goes here
%   Detailed explanation goes here

% 生成开边界条件
hin = obj.h0;
uin = obj.u0;
qxin = hin*uin;
% 读取边界节点编号
vert = read_edge_file(casename);
Nv = numel(vert); % 顶点个数
time = 0; % 开边界数据时间
Nt = numel(time); % 开边界文件时间步数
Nfield = 3; % 开边界物理场个数
f_Q = zeros(Nv, Nfield, Nt);
f_Q(:, 1) = hin;
f_Q(:, 2) = qxin;

% 生成开边界文件
obcfile = ndg_lib.phys.obc_file();
obcfile.make_obc_file([casename, '.nc'], time, vert, f_Q);
obcfile.set_file([casename, '.nc']);

end


function vert = read_edge_file(casename)
filename = [casename, '.edge'];
fp = fopen(filename);
% read total number of 
Nf = fscanf(fp, '%d', 1);
fgetl(fp); % pass the rest of first line

% read face to vertex list
data = fscanf(fp, '%d %d %d %d', [4, Nf]);
ind = data([2,3], :);
vert = unique( ind(:) );

fclose(fp);
end% func

