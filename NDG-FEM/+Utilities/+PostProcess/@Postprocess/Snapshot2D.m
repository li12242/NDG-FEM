function ph = Snapshot2D(obj, varname, stime, fileID ,varargin)
% Usages:
% 
x   = obj.NcFile(fileID).GetVarData('x');
y   = obj.NcFile(fileID).GetVarData('y');
[np, ne] = size(x);
var      = obj.GetVarData(varname, stime, fileID);
%var(var<1e-2) = nan;
vertex   = [x(:), y(:), var(:)];
bclist   = obj.StdCell.verlist';
EToV     = ones(ne, 1)*bclist;
EToV     = EToV + (np*(0:ne-1))'*ones(size(bclist));

switch varargin{1}
    case 'value' % ���ݺ���ֵ������ɫ
        % varargin ���Ӳ���??2.��ɫ1 3.��ɫ2 4.ֵ��Χ 5.����patch�������
        
        % colormap
        np = 40;
        color1 = varargin{2}; % .e.g, red - [1, 0, 0]
        color2 = varargin{3};
        map = [linspace(color1(1), color2(1), np); 
            linspace(color1(2), color2(2), np);
            linspace(color1(3), color2(3), np)]';
        colormap(map)
        
        % reinitialize clim
        set(gca, 'Clim', varargin{4});
        
        ph = patch(...
            'Vertices', vertex, ...
            'Faces', EToV, ...
            'FaceColor', 'interp', ...
            'FaceVertexCData', var(:),...
            varargin{5:end});  
        
    otherwise % Ĭ��
        ph = patch('Vertices', vertex, 'Faces', EToV, ...
            'FaceColor', [0.8, 0.9, 1]);
end
end% func
