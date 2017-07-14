        function nodelist = GetFaceListToNodeList(obj)
            % return the node list of face node
            % size [nFaceNode x 1]
            nodelist = zeros(obj.nFaceNode, 1);
            nodelist(1) = 1;
            nodelist(2) = obj.nOrder+1;
            nodelist = int16(nodelist);
        end %function