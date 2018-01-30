function [NcMax, Ncv, VToK, VToM] = assembleVertexCellConnect( obj )

Nv = obj.Nv; % total number of vertex
Ncv = zeros(Nv, 1); % number of cells connecting to each vertex
for m = 1:obj.Nmesh
    for k = 1:obj.meshUnion(m).K
        v = obj.meshUnion(m).EToV(:, k); % get the vertex index
        Ncv(v) = Ncv(v) + 1;
    end
end
% get the max number of cells connecting to each vertex
NcMax = max(Ncv);

VToK = zeros(NcMax, Nv);
VToM = zeros(NcMax, Nv);
Ncv = zeros(Nv, 1);
for m = 1:obj.Nmesh
    for k = 1:obj.meshUnion(m).K
        v = obj.meshUnion.EToV(:, k); % get the vertex index
        ind = Ncv(v)+1 + (v-1)*NcMax;
        VToK(ind) = k;
        VToM(ind) = m;
        Ncv(v) = Ncv(v) + 1;
    end
end

VToW = zeros(NcMax, Nv);
for n = 1:Nv
    w = zeros(Ncv(n), 1);
    for m = 1:Ncv(n)
        cellId = VToK(m, n);
        meshId = VToM(m, n);
        xc = obj.meshUnion(meshId).xc(cellId);
        yc = obj.meshUnion(meshId).yc(cellId);
        w(m) = 1./( (obj.meshUnion(meshId).vx(n) - xc).^2 ...
            + (obj.meshUnion(meshId).vy(n) - yc).^2 );
    end
    VToW( 1:Ncv(n), n ) = w./sum(w);
end
end% func