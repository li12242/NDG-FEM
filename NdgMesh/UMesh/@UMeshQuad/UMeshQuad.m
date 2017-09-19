classdef UMeshQuad < UMeshUnion2d
    %QUAD_MESH 四边形网格类
    %   根据输入参数生成对应四边形网格对象，共有三种生成方法
    %   1. 'file' 通过网格文件生成。
    %       输入参数为文件名（不包括后缀），所有文件包括单元节点文件，节点坐标文件，
    %       边界节点文件。文件格式见
    %       'ndg_lib/mesh/mesh2d/private/read_from_file.m'
    %       
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
    %       调用四边形均匀网格生成程序，输入参数包括：
    %           xlim，ylim - x，y 坐标范围；
    %           Mx，My - x 与 y 方向单元个数；
    %           facetype - 底部、上部、左侧和右侧边界类型；
    %
    
    methods(Static)
        [Nv, vx, vy, vz, K, EToV, EToR, EToBS] ...
            = uniform_mesh(xlim, ylim, zlim, Mx, My, Mz, facetype)
    end
    
    methods
        function obj = UMeshQuad(cell, varargin)
            % check input element
            if( ne(cell.type, ndg_lib.std_cell_type.Quad) )
                error(['Input cell type ', cell.type, ...
                    'is not quadrilateral!'])
            end
            
            obj = obj@UMeshUnion2d(cell, varargin{:});
                        
        end% func
        
        obj = refine(obj, multi_rate); % 加密网格
    end% methods
    
end

