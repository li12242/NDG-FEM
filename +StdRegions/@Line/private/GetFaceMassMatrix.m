function [ Mes ] = GetFaceMassMatrix( obj )
         Mes = zeros(obj.nNode,obj.nFaceNode);
        for ib = 1:obj.nFace
                % iFaceList: the No. of ibth boundary local face node list. 
                iFaceList = obj.GetFaceListAtFace(ib);
                [~,iFaceNodeList] = obj.GetNodeListAtFace(ib);
                Mes(iFaceNodeList,iFaceList) = 1;
        end
end