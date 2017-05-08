function [Nedge, Nnode, kM, kP, fM, fP, ftype, idM, idP, fpM, fpP, fscal, ...
    fnxM, fnyM, fnzM] = edge_connect(obj, EToV, EToE, EToF, EToBS)
%MESH_EDGE_CONNECT Summary of this function goes here
%   Detailed explanation goes here

% count the unique edges in mesh
[Nedge, kM, kP, fM, fP, ftype] = face_connect(obj, EToV, EToE, EToF, EToBS);

% collect the local and adjacent node index of each edge nodes
[Nnode, idM, idP, fpM, fpP, fscal, fnxM, fnyM] = ...
    node_connect(obj, Nedge, kM, kP, fM, fP);

fnzM = zeros(size(fnxM));
end% func


function [Nedge, kM, kP, fM, fP, ftype] = ...
    face_connect(obj, EToV, EToE, EToF, EToBS)

% count the unique edges for triangle
ind = zeros(obj.cell.Nface, obj.K);
for f = 1:obj.cell.Nface
    
    % find two vertex index of each edge
    v1 = EToV(obj.cell.FToV(1,f), :);
    v2 = EToV(obj.cell.FToV(2,f), :);
    % calculate the indicator for each edge
    ind(f, :) = min(v1, v2)*obj.Nv + max(v1, v2);
end
% calculate the unique edges
[~, id, ~] = unique(ind); % find the global index

Nedge = numel(id);
kM = fix( (id-1)./obj.cell.Nface )+1; % trans global index to element index
fM = rem(id-1, obj.cell.Nface)+1; % trans global index to local face index
kP = EToE(id); % adjacent element index
fP = EToF(id); % adjacent face index
ftype = EToBS(id); % face type
end% func

function [Nnode, idM, idP, fpM, fpP, fscal, fnxM, fnyM] = ...
    node_connect(obj, Nedge, kM, kP, fM, fP)

% calculate the total nodes
Nfp = obj.cell.N + 1; % nodes on each edge (line)
Nnode = Nedge*Nfp;

idM = zeros(Nnode, 1); idP = zeros(Nnode, 1);
fpM = zeros(Nnode, 1); fpP = zeros(Nnode, 1);
fscal = zeros(Nnode, 1);
fnxM = zeros(Nnode, 1); fnyM = zeros(Nnode, 1);

Np = obj.cell.Np; % number in each cell
Fmask = obj.cell.Fmask; % local index (elemental) of each face node 
faceIndexStart = zeros(obj.cell.Nface, 1); % start index of each face node
for f = 2:obj.cell.Nface
    faceIndexStart(f) = faceIndexStart(f-1) + obj.cell.Nfp(f-1);
end
sk = 1;
list = 1:Nfp;
for f = 1:Nedge
    k1 = kM(f); k2 = kP(f);
    f1 = fM(f); f2 = fP(f);

    ind2 = (k2-1)*Np + Fmask(:, f2);
    for n = 1:Nfp
        idM(sk) = (k1-1)*Np + Fmask(n, f1);
        fpM(sk) = faceIndexStart(f1)+n;
        xpM = obj.x(idM(sk));
        ypM = obj.y(idM(sk));
        
        xP = obj.x( ind2 );
        yP = obj.y( ind2 );
        d12 = (xpM - xP).^2 + (ypM - yP).^2;
        m = (d12 < 1e-10);
        idP(sk) = (k2-1)*Np + Fmask(m, f2);
        fpP(sk) = faceIndexStart(f2)+list(m);
        
        fscal(sk) = obj.Js(f1, k1)./obj.J( idM(sk)  );
        fnxM(sk) = obj.nx(f1, k1);
        fnyM(sk) = obj.ny(f1, k1);
        sk = sk+1;
    end
end% for
end
