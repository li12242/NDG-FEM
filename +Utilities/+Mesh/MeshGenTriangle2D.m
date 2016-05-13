function [VX,VY,EToV] = MeshGenTriangle2D(ne,start_coor,end_coor,flag)
% Mesh generater with triangle element
%
% DESCRIPTION
%   domain [start_coor, end_coor]x[start_coor, end_coor];
%   ne - No. of element on each edge
%
%   index
%       x coordinate 
% 
%       1   2   3   4   5   ....            n  n+1
%       -----------------------------------------
%       |   |   |   |   |   |   |   |   |   |   |
% (n+2) ----------------------------------------- 2(n+1)
%
% INPUT
%    ne      = No. of elements on each edge
%    start_coor = start coordinate
%    end_coor   = end coordinate
%    flag   = partination type [flag=0, "\", flag=1, "/"]
%
% OUTPUT
%    X = x coordinate
%    Y = y coordinate
%    E = index of vertices in element
%
% EXAMPLE USAGE
%    [X,Y,E] = MeshGenTriangle2D(21, 0, 1, 1)
%
% Author(s)
%    li12242 Tianjin University
%
% Reversion
%    v1.0 2014-12-8
%========================================================================== 

% info of point
x = linspace(start_coor, end_coor, ne+1); y = linspace(end_coor, start_coor, ne+1);
x = repmat(x, ne+1, 1); y = repmat(y, ne+1, 1);
x = x';
VX = x(:); VY = y(:);

EToV = zeros(2*ne^2,3);
for i = 1:ne % each row
    upLayerNum = (ne+1)*(i-1)+1:(ne+1)*(i-1)+ne+1;
    downLayerNum = upLayerNum + (ne+1);

    if flag     %
        EToV(2*ne*(i-1)+1:2*ne*(i-1)+ne,:) = ...
            [downLayerNum(1:ne)', upLayerNum(2:ne+1)', upLayerNum(1:ne)'];
        EToV(2*ne*(i-1)+ne+1:2*ne*i,:) = ...
            [downLayerNum(1:ne)',downLayerNum(2:ne+1)',upLayerNum(2:ne+1)'];
        
    else
        EToV(2*ne*(i-1)+1:2*ne*(i-1)+ne,:) = ...
            [downLayerNum(2:ne+1)', upLayerNum(2:ne+1)', upLayerNum(1:ne)'];
        EToV(2*ne*(i-1)+ne+1:2*ne*i,:) = ...
            [downLayerNum(1:ne)',downLayerNum(2:ne+1)',upLayerNum(1:ne)'];
        
    end

end
end