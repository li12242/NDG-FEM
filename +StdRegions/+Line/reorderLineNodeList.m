%% reorderLineNodeList


function facelist = reorderLineNodeList(nOrder, vOrder)
% 按照本单元(this)节点顺序，返回对方单元(next)节点编号
% 
% Input: 
%   nOrder  - order of polynomial
%   vOrder  - 节点顺序, nextVertice(vOrder) = thisVertice
% Output:
%   facelist    - 
% 
% 使用 thisElement 标准单元局部坐标作为直线顶点坐标
VerCoor = [-1; 1];

% nextElement 对应顶点坐标
nextVerCoor = zeros(size(VerCoor));
nextVerCoor(vOrder) = VerCoor(:);

% 获得 thisElement 标准单元节点坐标
thisNodeCoor = Node1D(nOrder+1);
% 将 thisElement 节点坐标 投影到 nextElement 单元
NodeCoor = (1-thisNodeCoor)/2.*nextVerCoor(1) + (1+thisNodeCoor)/2.*nextVerCoor(2);

% 根据节点坐标计算函数 $f = r$
nextf = NodeCoor;

% 在单元 thisElement 函数值 f 随节点顺序逐渐增大
% 单元 thisElement 节点与 nextElement 相互对应, 对应节点 f 值相同
% 根据 f 值大小可得本单元节点编号, 重排 f 获得本单元节点顺序
% thisf = nextf(facelist)
[thisf, facelist] = sort(nextf);
end% func

function Coor = Node1D(nNode)
Coor = linspace(-1, 1, nNode);
end% func