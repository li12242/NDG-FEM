function [ Nfp, nodelist ] = GetNodeListAtFace( obj, ind )
%GETNODELISTATFACE ��ȡ��ind�����Ͻڵ�ȫ�ֽڵ�����ڵ����
%   �ڵ㰴����ʱ������

order = obj.nOrder;
Nfp = order+1;
facelist = obj.GetFaceListAtFace(ind);
faceListToNodeList = obj.GetFaceListToNodeList();
nodelist = faceListToNodeList(facelist);

end