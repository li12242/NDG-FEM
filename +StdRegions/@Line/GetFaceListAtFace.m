    function facelist = GetFaceListAtFace(~, iface) 
            % return the face list at spicific face
            if (iface ==1)
                contour = 0; % the sum of boundary node num before iface
            else 
                contour =1;
            end
            % the number of nodes at the ibth boundary is 'nPerBoundaryNode(ib)'
            facelist = contour+1:contour+1;
        end %function getLocalFaceList