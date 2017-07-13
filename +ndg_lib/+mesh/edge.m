classdef edge
    %EDGE 网格内边界对象
    %   Detailed explanation goes here
    
    properties
        cell % 标准单元
        mesh % 网格对象
    end
    
    properties(SetAccess=protected)
        Nedge % 边界个数
        Nnode % 边界节点个数
        kM, kP % 相邻两侧单元编号
        fM, fP % 所在相邻两个单元面编号
        ftype@int8 % 面类型
        idM, idP % 边界上两侧节点编号
        fpM, fpP % 边界上两侧节点在面节点中局部编号
        fJs % 雅克比行列式
        fnxM, fnyM, fnzM % 法向量
    end
    
    methods
        function obj = edge(mesh)
            obj.cell = mesh.cell;
            obj.mesh = mesh;
            obj = obj.edge_connect(mesh.EToE, mesh.EToF, mesh.EToBS);
            obj = node_connect(obj, mesh, ...
                obj.Nedge, obj.kM, obj.kP, obj.fM, obj.fP);
        end
    end
    
    %% 私有函数
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

