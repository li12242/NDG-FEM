function [ fdr, fds ] = GradOrthogonalFun( r,s,order,ind )
%GRADORTHOGONALFUN ��׼��������������������ֵ

% ����ͶӰ������
[a,b] = rstoab(r,s);

% ���������ת��Ϊ�����ں�����ţ�i��j��
[i, j] = TransInd(order,ind);

% ����ԣ�r,s������ĵ���
[fdr, fds] = GradSimplex2DP(a,b,i,j);
end