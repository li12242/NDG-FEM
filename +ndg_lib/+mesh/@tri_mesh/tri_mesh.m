classdef tri_mesh < ndg_lib.mesh.mesh2d
    %TRI_MESH triangular mesh object
    %   根据输入参数生成对应四边形网格对象，共有三种生成方法
    %   1. 'file' 通过网格文件生成。
    %       输入参数为文件名（不包括后缀），所有文件包括单元节点文件，节点坐标文件，
    %       边界节点文件。文件格式见
    %       'ndg_lib/mesh/mesh2d/private/read_from_file.m'
    %
    %   2. 'variable' 通过变量生成网格对象。
    %       按照顺序给定各个参数：
    %           Nv - 顶点个数；
    %           vx，vy，vz - 顶点坐标；
    %           K - 单元个数；
    %           EToV - 单元顶点编号；
    %           EToR - 单元类型；
    %           EToBS - 单元边界类型；
    %   3. 'uniform' 生成均匀网格对象。
    %       调用三角形均匀网格生成程序，输入参数包括：
    %           xlim，ylim - x，y 坐标范围；
    %           Mx，My - x 与 y 方向单元个数；
    %           facetype - 底部、上部、左侧和右侧边界类型；
    %
    
    methods(Static)
        [Nv, vx, vy, vz, K, EToV, EToR, EToBS] ...
            = uniform_mesh(xlim, ylim, zlim, Mx, My, Mz, facetype)
    end
    
    methods
        function obj = tri_mesh(cell, varargin)
            switch varargin{1}
                case 'file'
                    [Nv, vx, vy, K, EToV, EToR, EToBS] ...
                        = read_from_file( varargin{2} );
                case 'variable'
                    var = varargin{2};
                    [Nv, vx, vy, K, EToV, EToR, EToBS] ...
                        = check_user_input(var);
                case 'uniform'
                    var = varargin{2};
                    xlim = var{1};
                    ylim = var{2};
                    Mx = var{3};
                    My = var{4};
                    facetype = var{5};
                    [Nv, vx, vy, K, EToV, EToR, EToBS] ...
                        = uniform_mesh(xlim, ylim, Mx, My, facetype);
                otherwise
                    
            end
            obj = obj@ndg_lib.mesh.mesh2d(cell, ...
                Nv, vx, vy, K, EToV, EToR, EToBS);
            
            % check input element
            if( ne(obj.cell.type, ndg_lib.std_cell_type.Tri) )
                error(['Input cell type ', cell.type, 'is not triangle!'])
            end
            
        end% func
        
        obj = refine(obj, multi_rate); % 加密网格
    end% methods
    
end

function [Nv, vx, vy, K, EToV, EToR, EToBS] = check_user_input(input)
% check the input variables for initlizing the mesh object.
Nv = input{1};
vx = input{2}; 
vy = input{3}; 
K  = input{4}; 
EToV = input{5};
EToR = int8(input{6}); 
EToBS = int8(input{7});

% check the # of the vertex
if ( numel(vx) ~= Nv ) || ( numel(vy) ~= Nv )
    error(['The length of input vertex coordinate "vx" or "vy" ', ...
        'is not equal to Nv']);
end% func

% check the # of the elements
if ( size(EToV, 2) ~= K )||( size(EToBS, 2) ~= K ) ||( numel(EToR) ~= K )
    error(['The length of input "EToV", "EToR" or "EToBS" ', ...
        'is not equal to K']);
end% func

if ( size(EToV, 1) ~= size(EToBS, 1) )
    error('The numbers of vertex in "EToV" and "EToBS" are not equal');
end% func

end% func
