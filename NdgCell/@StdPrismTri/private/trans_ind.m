function [ i,j ] = trans_ind( N, ind )
%TRANS_IND 将正交基函数编号转换为矩形单元内正交基函数编号（i,j）
%   三角形内正交函数对应矩形单元内基函数顺序为
%   i = 0, j = 0,1,2,...,N-1,N;
%   i = 1, j = 0,1,2....,N-1;
%   ...
%   i = N, J = 0;

% N - 最高阶数
% ind - 基函数序号
sk = 1;
for i = 0:N
    for j = 0:(N-i)
        if (abs(sk-ind)<10e-4)
            return;
        end
        sk = sk+1;
    end
end% for
end% func
