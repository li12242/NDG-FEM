function [ fval ] = OrthogonalFun( r,s,order,ind )
%ORTHOGONALFUN �����ı����������������ڽڵ㴦����ֵ
%   �ı���������������ͨ��һά��������չ������

% ת�����Ϊ��i,j����ʽ
[i,j] = TransInd(order,ind);
% ������������������ֵ
fval = Polylib.JacobiP(r(:), 0, 0, i).*Polylib.JacobiP(s(:), 0, 0, j);

end

