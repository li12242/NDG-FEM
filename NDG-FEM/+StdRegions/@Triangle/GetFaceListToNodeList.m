function [ nodelist ] = GetFaceListToNodeList( obj )
%GETFACELISTTONODELIST 返回边界节点的全局（单元）节点编号
%   节点按照逆时针顺序排序

Nfp = obj.nOrder+1;   % 每边节点个数
nodelist = zeros(obj.nFaceNode,1);

array = Nfp:-1:1;
% face 1, s = -1
nodelist(1:Nfp) = 1:Nfp;
% face 3, r = -1
temp = ones(Nfp, 1);
for i = 2:Nfp
    temp(i) = temp(i - 1) + array(i-1);
end
nodelist(3*Nfp:-1:2*Nfp+1) = temp;
% face 2, r+s = 0
temp2 = Nfp*ones(Nfp, 1);
temp2(2:end-1) = temp(3:end) -1; 
temp2(end) = temp(end);
nodelist(Nfp+1:2*Nfp) = temp2;

end