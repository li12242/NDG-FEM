%> @brief Edge information on unstructured meshes.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDG-FEM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef UEdge < handle
    
    properties(Hidden=true, SetAccess=protected)
    end
    
    properties(SetAccess=protected)
        %> pointer to edge cell object
        bcell
        %> number of edges
        Nedge
        %> number of nodes
        Nnode
        %> vertex index on each edge
        FToV
        %> adjacent cell index [kM; kP]
        FToE
        %> local face index of adjacent cells [fM; fP]
        FToF
        %> local node index of local cells
        FToN1
        %> local node index of adjacent cells
        FToN2
        %> boundary condition types of edges
        bcType
        %> determination of Jacobian matrixs
        J
        %> x component of the outward normal vector
        nx
        %> y component of the outward normal vector
        ny
        %> z component of the outward normal vector
        nz
    end
    
    methods
        function obj = UEdge(cell, K, EToV, Nv, vx, vy, vz, BCToV)
            obj.cell = cell;
            obj = obj.edge_connect(mesh.EToE, mesh.EToF, mesh.EToBS);
            obj = node_connect(obj, mesh, ...
                obj.Nedge, obj.kM, obj.kP, obj.fM, obj.fP);
        end
    end
    
    methods(Hidden, Access=protected)
        function obj = edge_connect(obj, EToE, EToF, EToBS)
            % count the unique edges
            ind = obj.mesh.Eind;
            [~, id, ~] = unique(ind); % 寻找不同面标记号
            
            % 赋值
            obj.Nedge = numel(id);
            obj.kM = fix( (id-1)./obj.cell.Nface )+1; % 不同所b标记在列号
            obj.fM = rem(id-1, obj.cell.Nface)+1; % 不同标记所在行号
            obj.kP = EToE(id); % adjacent element index
            obj.fP = EToF(id); % adjacent face index
            obj.ftype = int8(EToBS(id)); % face type
            
        end
        
        function obj = node_connect(obj, mesh, Nedge, kM, kP, fM, fP)
            
            nnode = 0;
            for f = 1:Nedge % 计算所有节点个数之和
                nnode = nnode + obj.cell.Nfp(fM(f));
            end
            
            idm = zeros(nnode, 1);
            idp = zeros(nnode, 1);
            fpm = zeros(nnode, 1); % 面节点在左单元面节点集合中局部序号
            fpp = zeros(nnode, 1); % 面节点在右单元面节点集合中局部序号
            fjs = zeros(nnode, 1);
            fnxm = zeros(nnode, 1);
            fnym = zeros(nnode, 1);
            fnzm = zeros(nnode, 1);
            
            Np = obj.cell.Np;
            Fmask = obj.cell.Fmask;
            faceIndexStart = zeros(obj.cell.Nface, 1); % 每个面起始面节点编号
            for f = 2:obj.cell.Nface
                faceIndexStart(f) = faceIndexStart(f-1) + obj.cell.Nfp(f-1);
            end
            
            sk = 1;
            for f = 1:Nedge
                k1 = kM(f); k2 = kP(f);
                f1 = fM(f); f2 = fP(f);
                Nfp = obj.cell.Nfp(f1);
                list = 1:Nfp;
                
                ind2 = (k2-1)*Np + Fmask(:, f2);
                for n = 1:Nfp
                    idm(sk) = (k1-1)*Np + Fmask(n, f1);
                    fpm(sk) = faceIndexStart(f1)+n;
                    xpM = mesh.x(idm(sk));
                    ypM = mesh.y(idm(sk));
                    zpM = mesh.z(idm(sk));
                    
                    xP = mesh.x( ind2 );
                    yP = mesh.y( ind2 );
                    zP = mesh.z( ind2 );
                    d12 = (xpM - xP).^2 + (ypM - yP).^2 + (zpM - zP).^2;
                    m = (d12 < 3e-16);
                    %try
                    idp(sk) = (k2-1)*Np + Fmask(m, f2);
                    %catch
                    %    keyboard
                    %end
                    fpp(sk) = faceIndexStart(f2)+list(m);
                    
                    fjs(sk) = mesh.Js(fpm(sk), k1);
                    fnxm(sk) = mesh.nx(fpm(sk), k1);
                    fnym(sk) = mesh.ny(fpm(sk), k1);
                    fnzm(sk) = mesh.ny(fpm(sk), k1);
                    sk = sk+1;
                end
            end% for
            
            % 赋值
            obj.Nnode = nnode;
            obj.idM = idm;
            obj.idP = idp;
            obj.fpM = fpm;
            obj.fpP = fpp;
            obj.fJs = fjs;
            obj.fnxM = fnxm;
            obj.fnyM = fnym;
            obj.fnzM = fnzm;
        end
        
        
    end
end

