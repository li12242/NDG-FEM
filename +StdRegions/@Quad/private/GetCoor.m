function [ r,s ] = GetCoor( order )
%GETCOOR Summary of this function goes here
%   Detailed explanation goes here

np = order+1;
[x,~] = Polylib.zwglj(np);

% ���Ȱ���x����ѭ����Ȼ����y����ѭ��
r = x*ones(1, np);
s = ones(np, 1)*x';

r = r(:); s = s(:);
end

