function [ i,j ] = TransInd( N,ind )
%TRANSIND 将正交基函数编号转换为编号（i,j）
%   四边形内正交函数对应顺序为
%   i = 0, j = 0,1,2,...,N-1,N;
%   i = 1, j = 0,1,2....,N-1,N;
%   ...
%   i = N, j = 0,1,2....,N-1,N;

% N - 最高阶数
% ind - 基函数序号
sk = 1;
for i = 0:N
    for j = 0:N
        if (abs(sk-ind)<10e-4)
            return;
        end
        sk = sk+1;
    end
end% for
end

