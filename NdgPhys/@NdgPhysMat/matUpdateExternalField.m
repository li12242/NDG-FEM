%> Update the external field
function matUpdateExternalField( obj, time, fphys )

% % update external value from other mesh
% for m1 = 1:obj.Nmesh
%     mesh = obj.meshUnion( m1 );
%     
%     for e = 1:numel(mesh.edgeUnion)
%         edge = mesh.edgeUnion(e);
%         m2 = edge.FToM;
%         IntMat = edge.IntM;
%         for n = 1:edge.M
%             k1 = edge.FToE(1, n);
%             k2 = edge.FToE(2, n);
%             p1 = edge.FToN1(:, n);
%             p2 = edge.FToN2(:, n);
%             
%             for fld = 1:obj.Nfield
%                 obj.fext{m1}(p1, k1, fld) = IntMat * fphys{m2}(p2, k2, fld);
%             end
%         end
%     end
% end

end% func