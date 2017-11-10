function [Nvcmax, Nvc, VToK, VToM, VToW] = assembleVertexCellConnect( obj )

Nv = obj.Nv; % total number of vertex
Nvc = zeros(Nv, 1); % number of cells connecting to each vertex
for m = 1:obj.Nmesh
    for k = 1:obj.meshUnion(m).K
        v = obj.meshUnion(m).EToV(:, k); % get the vertex index
        Nvc(v) = Nvc(v) + 1;
    end
end
% get the max number of cells connecting to each vertex
Nvcmax = max(Nvc);

VToK = zeros(Nvcmax, Nv);
VToM = zeros(Nvcmax, Nv);
Nvc = zeros(Nv, 1);
for m = 1:obj.Nmesh
    for k = 1:obj.meshUnion(m).K
        v = obj.meshUnion.EToV(:, k); % get the vertex index
        ind = Nvc(v)+1 + (v-1)*Nvcmax;
        VToK(ind) = k;
        VToM(ind) = m;
        Nvc(v) = Nvc(v) + 1;
    end
end

VToW = zeros(Nvcmax, Nv);
for n = 1:Nv
    w = zeros(Nvc(n), 1);
    for m = 1:Nvc(n)
        cellId = VToK(m, n);
        meshId = VToM(m, n);
        xc = obj.meshUnion(meshId).xc(cellId);
        yc = obj.meshUnion(meshId).yc(cellId);
        w(m) = 1./( (obj.meshUnion(meshId).vx(n) - xc).^2 ...
            + (obj.meshUnion(meshId).vy(n) - yc).^2 );
    end
    VToW( 1:Nvc(n), n ) = w./sum(w);
end
end% func