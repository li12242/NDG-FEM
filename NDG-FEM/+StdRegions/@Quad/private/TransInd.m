function [ i,j ] = TransInd( N,ind )
%TRANSIND ���������������ת��Ϊ��ţ�i,j��
%   �ı���������������Ӧ˳��Ϊ
%   i = 0, j = 0,1,2,...,N-1,N;
%   i = 1, j = 0,1,2....,N-1,N;
%   ...
%   i = N, j = 0,1,2....,N-1,N;

% N - ��߽���
% ind - ���������
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

