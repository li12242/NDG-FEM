function [i,j] = TransInd(N, ind)
%TRANSIND ���������������ת��Ϊ���ε�Ԫ��������������ţ�i,j��
%   ������������������Ӧ���ε�Ԫ�ڻ�����˳��Ϊ
%   i = 0, j = 0,1,2,...,N-1,N;
%   i = 1, j = 0,1,2....,N-1;
%   ...
%   i = N, J = 0;

% N - ��߽���
% ind - ���������
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
