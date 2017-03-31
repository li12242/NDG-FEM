        function  [n, nodelist] = GetNodeListAtFace(obj, iface)
            
            facelist = obj.GetFaceListAtFace(iface);
            faceListToNodelist = GetFaceListToNodeList(obj);
            nodelist = faceListToNodelist(facelist);
            n =1;
        end% func