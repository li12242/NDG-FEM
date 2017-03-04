classdef Triangle < StdRegions.BaseElement
    %TRIANGLE ��ά��׼��Ԫ
    %   ��׼�����ε�Ԫ�����ýڵ��������nodal basis������Ӧ�ڵ�ΪLGL�ڵ�ֲ�
    %   ��Hesthaven and Warburton, 2008��
    
    properties
        nFaceNode % �߽�ڵ����
        r       % ��׼��Ԫ�ڵ�x����
        s       % ��׼��Ԫ�ڵ�y����
        V       % Vandermonde����
        invV    % Vandermonde�����
        M       % ��������
        Dr      % �ڵ��������r��
        Ds      % �ڵ��������s��
        Drw     % 
        Dsw
        Mes     % �߽�ڵ������������������
        LIFT    % �߽�ͨ��ת������
        Fmask   % �߽��Ͻڵ���
    end
    
    methods
        % ==========================
        % ���캯��
        % ==========================
        function obj = Triangle(order)
            dim = 2; vertice = 3; face = 3;
            % �̳����԰����� 
            % nDim, nVertice, nOrder, nFace, nNode 
            obj = obj@StdRegions.BaseElement(dim,vertice,order,face);
            % �̳������趨
            obj.sName = 'Triangle';
            obj.nNode = (order+1)*(order+2)/2;
            obj.nFaceNode = (order+1)*face;
            
            % ��ȡ�ڵ㼯��
            [obj.r, obj.s] = GetCoor(order);
            % �߽��Ͻڵ���
            
            % Vandermonde����
            obj.V = obj.GetVandMatrix(order,obj.r,obj.s);
            obj.invV = inv(obj.V);
            % �������� $M = (V*V^T)^{-1} = invV^T*invV$
            obj.M = obj.invV'*obj.invV;
            % ��������
            [obj.Dr, obj.Ds] = GetGrandMatrix( order,obj.r,obj.s,obj.V );
            
            % �߽���ֻ�������������
            obj.Mes = GetFaceMassMatrix( obj );
            
            % LIFT���� $M^{-1} \cdot Mes$
            obj.LIFT = (obj.V*(obj.V)')*obj.Mes;
            
            % Fmask ���д洢��ÿ�д洢ith�߽ڵ���
            nodelist = obj.GetFaceListToNodeList();
            obj.Fmask = reshape(nodelist, order+1, face);
        end% func
        
        % ==========================
        % ����Vandermonde����
        % ==========================
        V = GetVandMatrix(obj, order, r, s);
        
        % ==========================
        % ��ȡ��ind�����Ͻڵ�ȫ�ֽڵ�����ڵ����
        % ==========================
        [ Nfp, nodelist ] = GetNodeListAtFace( obj, ind );
        [ facelist ] = GetFaceListAtFace( obj,ind );
        
        [ nodelist ] = GetFaceListToNodeList( obj );
        [ nodelist ] = GetVertexNodeList( obj );
    end
    
end

