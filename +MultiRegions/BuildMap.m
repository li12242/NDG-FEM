function [obj] = BuildMap(obj, VX, VY, EToV, EToE, EToF)

NODETOL = 1e-10;
nFaceNode = obj.Shape.nFaceNode; % total number on all edges
vmapM = zeros(nFaceNode, obj.nElement);
vmapP = zeros(nFaceNode, obj.nElement);
for k1 = 1:obj.nElement
    vmapM(:,k1) = obj.Shape.getFaceListToNodeList + (k1-1)*obj.Shape.nNode;
end

one = ones(1, obj.Shape.nFaceNode./obj.Shape.nFace);
for k1 = 1:obj.nElement
    for f1 = 1:obj.Shape.nFace
        % find neighbor
        k2 = EToE(k1,f1); f2 = EToF(k1,f1);

        flist1 = obj.Shape.getFaceListAtFace(f1);
        flist2 = obj.Shape.getFaceListAtFace(f2);

        % reference length of edge
        v1 = EToV(k1,f1); v2 = EToV(k1, 1+mod(f1, obj.Shape.nFace));
        refd = sqrt( (VX(v1)-VX(v2))^2 + (VY(v1)-VY(v2))^2 );

        % find find volume node numbers of left and right nodes 
        vidM = vmapM(flist1, k1); vidP = vmapM(flist2, k2);
        x1 = obj.x(vidM); y1 = obj.y(vidM); x2 = obj.x(vidP); y2 = obj.y(vidP);
        x1 = x1*one;  y1 = y1*one;  x2 = x2*one;  y2 = y2*one;

        % Compute distance matrix
        D = (x1 -x2').^2 + (y1-y2').^2;
        [idM, idP] = find(sqrt(abs(D))<NODETOL*refd);
        vmapP(flist1(idM), k1) = vidP(idP);
    end
end
% Assignment
obj.vmapM = vmapM;
obj.vmapP = vmapP;
end% func
