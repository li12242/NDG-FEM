function [Np,r,s,t] = node_coor_func(obj, nOrder)
% Compute (x,y) nodes in equilateral triangle for polynomial of order
alpopt = [0.0000 0.0000 1.4152 0.1001 0.2751 0.9800 1.0999 ...
          1.2832 1.3648 1.4773 1.4959 1.5743 1.5770 1.6223 1.6258];  
% Set optimized parameter, alpha, depending on order N
if (nOrder<16)
    alpha = alpopt(nOrder);
else
    alpha = 5/3;
end
% total number of nodes
Np = (nOrder+1)*(nOrder+2)/2;
% Create equidistributed nodes on equilateral triangle
L1 = zeros(Np,1); L3 = zeros(Np,1); %L2 = zeros(Np,1); 
sk = 1;
for n=1:nOrder+1
    for m=1:nOrder+2-n
        L1(sk) = (n-1)/nOrder; L3(sk) = (m-1)/nOrder;
        sk = sk+1;
    end
end
L2 = 1.0-L1-L3;
x = -L2+L3; y = (-L2-L3+2*L1)/sqrt(3.0);
% Compute blending function at each node for each edge
blend1 = 4*L2.*L3; blend2 = 4*L1.*L3; blend3 = 4*L1.*L2;
% Amount of warp for each node, for each edge
warpf1 = warpfactor(nOrder,L3-L2); 
warpf2 = warpfactor(nOrder,L1-L3);
warpf3 = warpfactor(nOrder,L2-L1);
% Combine blend & warp
warp1 = blend1.*warpf1.*(1 + (alpha*L1).^2);
warp2 = blend2.*warpf2.*(1 + (alpha*L2).^2);
warp3 = blend3.*warpf3.*(1 + (alpha*L3).^2);
% Accumulate deformations associated with each edge
x = x + 1*warp1 + cos(2*pi/3)*warp2 + cos(4*pi/3)*warp3;
y = y + 0*warp1 + sin(2*pi/3)*warp2 + sin(4*pi/3)*warp3;
[r,s] = xytors(x,y);

t = zeros(size(r));
end

function warp = warpfactor(N, rout)
% Compute scaled warp function at order N based on rout interpolation nodes

% Compute LGL and equidistant node distribution
[LGLr,~] = zwglj(N+1); req = linspace(-1,1,N+1)';

% Compute V based on req
Veq = VandMatrix(N, req);
% Evaluate Lagrange polynomial at rout
Nr = length(rout); %Lmat=zeros(N+1,Nr); 
Pmat = zeros(N+1,Nr);
for i=1:N+1
  Pmat(i,:) = JacobiP(rout, 0, 0, i-1);
end
Lmat = Veq'\Pmat;
% Compute warp factor
warp = Lmat'*(LGLr - req);
% Scale factor
zerof = (abs(rout)<1.0-1.0e-10); sf = 1.0 - (zerof.*rout).^2;
warp = warp./sf;
end

function V = VandMatrix(N, r)
V = zeros(numel(r), N+1);
for j=0:N
    % P_{j-1}(r_i)$
    V(:,j+1) = JacobiP(r(:), 0, 0, j);
end% for
end% func

function [r,s] = xytors(x,y)
% From (x,y) in equilateral triangle to 
% (r,s) coordinates in standard triangle
L1 = (sqrt(3.0)*y+1.0)/3.0;
L2 = (-3.0*x - sqrt(3.0)*y + 2.0)/6.0;
L3 = ( 3.0*x - sqrt(3.0)*y + 2.0)/6.0;
r = -L2 + L3 - L1; s = -L2 - L3 + L1;
end