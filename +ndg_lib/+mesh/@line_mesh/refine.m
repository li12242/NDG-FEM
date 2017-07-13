function obj = refine( obj, refine_level )
%REFIND Summary of this function goes here
%   Detailed explanation goes here

for m = 1:refine_level
    
    [newNv, newVX] = new_vertex(obj);
    [newK, newEToV, newEToBS, newEToR] = new_elemet_info(obj);
    
    obj = ndg_lib.mesh.line_mesh(obj.cell, ...
        newNv, newVX, newK, newEToV, newEToR, newEToBS);
end

end

function [newNv, newVX] = new_vertex(obj)
newVX = obj.cell_mean(obj.x);
newVX = [obj.vx; newVX(:)];
newNv = obj.Nv + obj.K;
end

function [newK, newEToV, newEToBS, newEToR] = new_elemet_info(obj)
newK = obj.K*2;
newEToV = zeros(obj.cell.Nv, newK);
newEToBS = zeros(obj.cell.Nv, newK);
newEToR = int8(ones(newK, 1));

inner = ndg_lib.bc_type.Inner;
for k = 1:obj.K
    ksid = (k-1)*2;
    v = obj.EToV(:, k); % µ¥Ôª¶¥µã±àºÅ
    bc = obj.EToBS(:, k);
    
    newEToV(:, ksid+1) = [v(1), k+obj.Nv];
    newEToV(:, ksid+2) = [k+obj.Nv, v(2)];
    newEToBS(:, ksid+1) = [bc(1), inner];
    newEToBS(:, ksid+1) = [inner, bc(2)];
    newEToR(:, ksid+1) = obj.EToR(k);
    newEToR(:, ksid+2) = obj.EToR(k);
end
end% func
