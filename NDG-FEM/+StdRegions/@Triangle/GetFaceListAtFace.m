function [ facelist ] = GetFaceListAtFace( obj,ind )
%GETFACELISTATFACE ���ص�ind�����Ͻڵ�ֲ����
%   �ֲ����˳��ֻ���Ǳ߽��Ͻڵ㣬

Nfp = obj.nOrder + 1; % ÿ�߽ڵ����
contour = Nfp*(ind-1); % ǰn���߽ڵ�����
facelist = (contour+1):(contour+Nfp); % ��ind���߽ڵ���

end

