function dryEleFlag = SWE_DryEle2d(mesh, h, minDepth)
%% Determine the dry elements
% The dry elements are determined with the averaged water depth $\bar{h}_i$, 
% when the mean depth is less than $\bar{h}_i < h_{Dry}$, than the elemnt
% is determined as dry elements.

% Input : 
%   mesh     - mesh object
%   h        - water depth
%   minDepth - threadhold of water depth
% Output:
%   dryEleFlag - flag (1) for dry elements
% 

%% Evaluate the water volume
% The water volume is obtained as $V = \sum_{j=1}^{Np} J_j w_j h_j$ where 
% $h_j$ is the water depth on Gauss quadrature nodes, $w_j$ is the 
% quadrature weights. The quadrature weight $w_j$ can obtain from the mass
% matrix, where $w_j = \sum_{i=1}^{Np} M_{ij} = \sum_{i=1}^{Np} 
% \int_{\Omega_{std}} l_i l_j dA = \int_{\Omega_{std}} l_j dA$. Then the
% water volume is same as $V = \int_{\Omega_{std}} J \cdot h_h(x) dA
% = \int_{\Omega_{std}} J \cdot \sum_{j=1}^{Np} h_j l_j dA$
% 
M = mesh.Shape.M; AVE = sum(M);
% water volume
vol = AVE*(h.*mesh.J);

%% Get the mean water depth and dry element index

area  = AVE*mesh.J;
hmean = vol./area;
dryEleFlag = (hmean < minDepth);

end% func