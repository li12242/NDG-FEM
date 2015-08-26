function facelist = reorderTriNodeList(nOrder, vOrder)
% 按照本单元(this)节点顺序，返回对方单元(next)节点编号
% nextElementVertice(vOrder) = thisElementVertice
nextVerCoor.r = [-1; 1; -1];
nextVerCoor.s = [-1; -1; 1];

VerCoor.r = zeros(size(nextVerCoor.r));
VerCoor.s = zeros(size(nextVerCoor.s));
VerCoor.r(:) = nextVerCoor.r(vOrder);
VerCoor.s(:) = nextVerCoor.s(vOrder);

% next element coordiante
[x,y] = Nodes2D(nOrder);
[r, s] = xytors(x,y);
NodeCoor.r = 0.5*(-(r+s)*VerCoor.r(1)+(1+r)*VerCoor.r(2)+(1+s)*VerCoor.r(3));
NodeCoor.s = 0.5*(-(r+s)*VerCoor.s(1)+(1+r)*VerCoor.s(2)+(1+s)*VerCoor.s(3));
% NodeCoor.r & s is along "this" element node order

M = ceil(1./abs(NodeCoor.r(2)- NodeCoor.r(1)))+100;
thisf = NodeCoor.r + M* NodeCoor.s;

[~,facelist] = sort(thisf);
[~,facelist] = sort(facelist);
end

function [x,y] = Nodes2D(N)
% total number of nodes
Np = (N+1)*(N+2)/2;
% Create equidistributed nodes on equilateral triangle
L1 = zeros(Np,1); L3 = zeros(Np,1); %L2 = zeros(Np,1); 
sk = 1;
for n=1:N+1
  for m=1:N+2-n
    L1(sk) = (n-1)/N; L3(sk) = (m-1)/N;
    sk = sk+1;
  end
end
L2 = 1.0-L1-L3;
x = -L2+L3; y = (-L2-L3+2*L1)/sqrt(3.0);
end

function [r,s] = xytors(x,y)
% From (x,y) in equilateral triangle to 
% (r,s) coordinates in standard triangle
L1 = (sqrt(3.0)*y+1.0)/3.0;
L2 = (-3.0*x - sqrt(3.0)*y + 2.0)/6.0;
L3 = ( 3.0*x - sqrt(3.0)*y + 2.0)/6.0;
r = -L2 + L3 - L1; s = -L2 - L3 + L1;
end