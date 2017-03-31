classdef Line < StdRegions.BaseElement
      properties
        nFaceNode % 边界节点个数
        r       % 标准单元节点x坐标
       VandMatrix      % Vandermonde矩阵
        invV    % Vandermonde逆矩阵
        M       % 质量矩阵
        invM
        Dr      % 节点基函数对r求导
        Mes     % 边界节点基函数积分质量矩阵
        LIFT    % 边界通量转换矩阵
        Fmask   % 边界上节点编号
      end
        %%
    methods
       % ==========================
        % 构造函数
        % ==========================
        function obj=Line(order)
            dim=1;vertice=2;face=2;
            obj=obj@StdRegions.BaseElement(dim,vertice,order,face);
            obj.nFaceNode=face;
            obj.nNode=order+1;
            obj.sName='Line';
            obj.r=GetCoor(order);
            obj.VandMatrix=obj.GetVandMatrix(order,obj.r);
            obj.invV=inv(obj.VandMatrix);
            obj.M = (obj.invV)'*(obj.invV);
            obj.invM = obj.VandMatrix*(obj.VandMatrix)';
            obj.Dr=obj.GetDeriMatrix(obj.r);
            obj.Mes=GetFaceMassMatrix( obj );
            obj.LIFT=((obj.VandMatrix)*(obj.VandMatrix)')*obj.Mes;
            obj.Fmask=[1;order+1];
        end
    end
end


