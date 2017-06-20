function [ vertid ] = get_obc_vert( casename )
%GET_OBC_VERT 获取边界顶点序号
%   Detailed explanation goes here

filename = [casename, '.edge'];
fp = fopen(filename);
% read total number of 
Nf = fscanf(fp, '%d', 1);
fgetl(fp); % pass the rest of first line

% read face to vertex list
data = fscanf(fp, '%d %d %d %d', [4, Nf]);
surfid = data(4, :);
ind = ( surfid > 3 );
vert = data([2,3], ind);
vertid = unique( vert(:) );

fclose(fp);
end

