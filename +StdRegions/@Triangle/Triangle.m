classdef Triangle < StdRegions.BaseElement
    %TRIANGLE 二维标准单元
    %   标准三角形单元，采用节点基函数（nodal basis），对应节点为LGL节点分布
    %   （Hesthaven and Warburton, 2008）
    
    properties
        r       % 标准单元节点x坐标
        s       % 标准单元节点y坐标
        V       % Vandermonde矩阵
        invV    % Vandermonde逆矩阵
        M       % 质量矩阵
        Dr      % 节点基函数对r求导
        Ds      % 节点基函数对s求导
        Drw     % 
        Dsw
        Mes     % 边界积分节点质量矩阵
        LIFT    % 边界通量转换矩阵
        Fmask   % 边界上节点编号
    end
    
    methods
        function obj = Triangle(order)
            dim = 2; vertice = 3; face = 3;
            % 继承属性包括： 
            % nDim, nVertice, nOrder, nFace, nNode 
            obj = obj@StdRegions.BaseElement(dim,vertice,order,face);
            % 继承属性设定
            obj.sName = 'Triangle';
            obj.nNode = (order+1)*(order+2)/2;
            
            % 获取节点集合
            [obj.r, obj.s] = GetCoor(order);
            % 边界上节点编号
            
            % Vandermonde矩阵
            obj.V = obj.GetVandMatrix(order,obj.r,obj.s);
            obj.invV = inv(obj.V);
            % 质量矩阵
            obj.M = obj.invV'*obj.invV;
            % 导数矩阵
            [obj.Dr, obj.Ds] = GetGrandMatrix( order,obj.r,obj.s,obj.V );
            
        end% func
        
        % 计算Vandermonde矩阵
        V = GetVandMatrix(obj, order, r, s);
    end
    
end

