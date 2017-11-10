function draw( obj, varargin )
%DRAW Summary of this function goes here
%   Detailed explanation goes here

switch nargin
    case 1 % 
        f = obj.f_Q(:,:,1);
    case 2 % 标量场作为第二个参数输入
        f = varargin{1};
end
if ( isempty(obj.draw_h) || ~isvalid(obj.draw_h))
    % 若图像未绘制或窗口被关闭
    EToV = ones(obj.mesh.K, 1)*obj.mesh.cell.Fmask(:)';
    EToV = EToV + ( obj.mesh.cell.Np*(0:obj.mesh.K-1) )'...
        *ones(1, obj.mesh.cell.Nfptotal);
    obj.draw_h = patch(...
        'Vertices', [obj.mesh.x(:), obj.mesh.y(:), f(:)], ...
        'Faces', EToV, ...
        'FaceColor', 'interp', ...
        'FaceVertexCData', f(:));
else % 若图像存在
    set(obj.draw_h, ...
        'Vertices', [obj.mesh.x(:), obj.mesh.y(:), f(:)],...
        'FaceVertexCData', f(:));
end
end

