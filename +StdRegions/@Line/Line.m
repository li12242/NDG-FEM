classdef Line < StdRegions.BaseElement
      properties
        nFaceNode % �߽�ڵ����
        r       % ��׼��Ԫ�ڵ�x����
       VandMatrix      % Vandermonde����
        invV    % Vandermonde�����
        M       % ��������
        invM
        Dr      % �ڵ��������r��
        Mes     % �߽�ڵ������������������
        LIFT    % �߽�ͨ��ת������
        Fmask   % �߽��Ͻڵ���
      end
        %%
    methods
       % ==========================
        % ���캯��
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


