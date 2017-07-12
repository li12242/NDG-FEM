function refine_mesh( obj, multi_rate )
%REFINE_MESH Summary of this function goes here
%   Detailed explanation goes here

for m = 1:multi_rate
    xc = obj.mesh.cell_mean(obj.mesh.x); % 计算单元形心
    yc = obj.mesh.cell_mean(obj.mesh.y);
    xf = obj.mesh.face_mean(obj.mesh.x); % 计算单元边中点
    yf = obj.mesh.face_mean(obj.mesh.y);
    xe = zeros(obj.mesh.Nedge, 1);
    ye = zeros(obj.mesh.Nedge, 1);
    for n = 1:obj.mesh.Nedge
        xe(n) = xf(obj.mesh.fM(n), obj.mesh.kM(n));
        ye(n) = yf(obj.mesh.fM(n), obj.mesh.kM(n));
    end
    newNv = obj.mesh.Nv + obj.mesh.Nedge + obj.mesh.K;
    vx = [obj.mesh.vx; xe; xc'];
    vy = [obj.mesh.vy; ye; yc'];
    
    FToE = zeros(obj.mesh.cell.Nface, obj.mesh.K); % 每个单元周围面编号
    EToF = zeros(obj.mesh.cell.Nface, obj.mesh.K); % 每个单元周围面顺序
    tid = ones(obj.mesh.K, 1); % 每个单元计数器
    for n = 1:obj.mesh.Nedge
        kM = obj.mesh.kM(n);
        FToE(tid(kM), kM) = n + obj.mesh.Nv;
        EToF(tid(kM), kM) = obj.mesh.fM(n);
        tid(kM) = tid(kM) + 1;
        kP = obj.mesh.kP(n);
        if kM ~= kP
            FToE(tid(kP), kP) = n + obj.mesh.Nv;
            EToF(tid(kP), kP) = obj.mesh.fP(n);
            tid(kP) = tid(kP) + 1;
        end
    end
    CNToE = (1:obj.mesh.K) + obj.mesh.Nv + obj.mesh.Nedge;
    newK = 4*obj.mesh.K;
    newEToV = zeros(obj.mesh.cell.Nv, newK);
    newEToBS = zeros(obj.mesh.cell.Nv, newK);
    for k = 1:obj.mesh.K
        ksid = (k-1)*4 + 1;
        [~, ind] = sort(EToF(:, k));
        v = obj.mesh.EToV(:, k);
        bc = obj.mesh.EToBS(:, k);
        inner = ndg_lib.bc_type.Inner;
        
        newEToV(:, ksid  ) = [v(1), FToE(ind(1),k), CNToE(k), FToE(ind(4), k)]';
        newEToV(:, ksid+1) = [v(2), FToE(ind(2),k), CNToE(k), FToE(ind(1), k)]';
        newEToV(:, ksid+2) = [v(3), FToE(ind(3),k), CNToE(k), FToE(ind(2), k)]';
        newEToV(:, ksid+3) = [v(4), FToE(ind(4),k), CNToE(k), FToE(ind(3), k)]';
        newEToBS(:, ksid  ) = [bc(1), inner, inner, bc(4)];
        newEToBS(:, ksid+1) = [bc(2), inner, inner, bc(1)];
        newEToBS(:, ksid+2) = [bc(3), inner, inner, bc(2)];
        newEToBS(:, ksid+3) = [bc(4), inner, inner, bc(3)];
    end
    newEToR = int8(ones(newK, 1)).*ndg_lib.mesh_type.Normal;
    
    newMesh = ndg_test.mesh_test.mesh2d_fullquad(obj.mesh.cell, ...
        newNv, vx, vy, ...
        newK, newEToV, newEToR, newEToBS);
    obj.mesh = newMesh;
end

end

