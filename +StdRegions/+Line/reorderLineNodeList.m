function facelist = reorderLineNodeList(nOrder, vOrder)
% 按照本单元(this)节点顺序，返回对方单元(next)节点编号
% nextVertice(vOrder) = thisVertice
nextVerCoor = [-1; 1];
thisVerCoor = zeros(size(nextVerCoor));
thisVerCoor(:) = nextVerCoor(vOrder);

% [nextNodeCoor,~] = Polylib.zwglj(nOrder+1);
nextNodeCoor = Node1D(nOrder+1);
NodeCoor = (1-nextNodeCoor)/2.*thisVerCoor(1) + (1+nextNodeCoor)/2.*thisVerCoor(2);

f = NodeCoor;
[~,facelist] = sort(f);
[~,facelist] = sort(facelist);
end

function Coor = Node1D(nNode)
Coor = linspace(-1, 1, nNode);
end