%% SL2 slope limiter
% Slope limiter proposed by Anastasion and Chan (1997). 
% For more details, refer to Khan and Lai (2014).
function hlim = SLLoc2(mesh, h, beta)
% Input:
%   beta - Parameter ranges from 1 to 2 result in the minmod limiter and 
%           Superbee limiter, respectively.
% Output:
%   hlim - limited values
% 
%% Usages
%   shape        = StdRegions.Triangle(1);
%   [VX,VY,EToV] = Utilities.Mesh.MeshGenTriangle2D(4,4,-1,1,-1,1,0);
%   mesh         = MultiRegions.RegionTri(shape, EToV, VX, VY);
%   h            = (mesh.x + mesh.y);
%   beat         = 2; % superbee limiter
%   hlim         = Utilities.Limiter.Limiter2D.SL2(mesh, h, beat);
% 

% shape  = mesh.Shape;
% AVE    = sum(shape.M);
% Area   = AVE*mesh.J; % area of each element
% % AVE    = AVE./Area;
% %% Get the cell average
% hmean  = AVE*(mesh.J.*h)./Area;
% 
% %% Get the unlimited gradient
% % 
% % # Compute the value at the boundary
% % # Compute the unlimited gradient by Green¡¯s theorem
% % 
% % $\begin{array}{l}
% % \frac{\partial u}{\partial x} = \frac{1}{\Omega_0} \oint u \rm{dy} \cr
% % \frac{\partial u}{\partial y} = -\frac{1}{\Omega_0} \oint u \rm{dx}
% % \end{array}$
% % 
% 
% % 1.Compute the value at the boundary
% hb     = zeros(shape.nFace, mesh.nElement);
% % find the boundary value by the vertex value
% hv     = h(shape.getVertexNodeList,:); % vertex values
% for f1 = 1:shape.nFace
%     v1 = f1; v2 = mod(f1, shape.nFace)+1;
%     hb(f1, :) = (hv(v1, :) + hv(v2, :))/2;
% end% for
% 
% % 2.Compute the unlimited gradient by Green¡¯s theorem
% % vertex coordinate
% xv     = mesh.x(shape.getVertexNodeList,:);
% yv     = mesh.y(shape.getVertexNodeList,:);
% % initialize gradient
% hpx    = zeros(1, mesh.nElement); % initialize x gradient 
% hpy    = zeros(1, mesh.nElement); % initialize y gradient 
% % 
% for f1 =1:shape.nFace
%     v1  = f1; v2 = mod(f1, shape.nFace)+1;
%     dx  = (xv(v2,:) - xv(v1,:));
%     dy  = (yv(v2,:) - yv(v1,:));
%     hpx = hpx + hb(f1, :).*dy;
%     hpy = hpy - hb(f1, :).*dx;
% end% for
% hpx    = hpx./Area;
% hpy    = hpy./Area;
% 
% %% Compute the limited gradient
% hmax   = ones(shape.nNode, 1)*max([hmean(mesh.EToE'); hmean]); 
% hmin   = ones(shape.nNode, 1)*min([hmean(mesh.EToE'); hmean]); 
% h0     = ones(shape.nNode, 1)*hmean;
% 
% % 
% r      = ones(size(h));
% flag   = h>hmax;
% r(flag)= (hmax(flag) - h0(flag))./(h(flag) - h0(flag));
% flag   = h<hmin;
% r(flag)= (hmin(flag) - h0(flag))./(h(flag) - h0(flag));
% 
% varphi = max(min(beta*r, 1), min(r, beta));
% phi    = min(varphi);  % limited coefficient in each element
% %% Reconstruction
% % Reconstruct the nodes value
% xc     = (AVE*(mesh.J.*mesh.x))./Area; % central coordinate
% yc     = (AVE*(mesh.J.*mesh.y))./Area; % central coordinate
% 
% phpx    = hpx.*phi;
% phpy    = hpy.*phi;
% 
% dx      = mesh.x - ones(shape.nNode, 1)*xc;
% dy      = mesh.y - ones(shape.nNode, 1)*yc;
% 
% hlim    = h0 + dx.*(ones(shape.nNode,1)*phpx) ...
%     + dy.*(ones(shape.nNode,1)*phpy);

%% mex 
hlim = Utilities.Limiter.Limiter2D.SLLoc2d_Mex(h, mesh.J, mesh.Shape.M,...
    mesh.Shape.Fmask, mesh.EToE, mesh.x, mesh.y, beta);

end% func
