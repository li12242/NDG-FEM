classdef mesh2d_fullquad < ndg_lib.mesh.mesh2d
    %MESH2D_FULLQUAD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        invM % inverse mass matrix of each element
        Jq, % 插值节点处雅克比行列式
        rxq, ryq, rzq
        sxq, syq, szq
        txq, tyq, tzq
        
        eidMq, eidPq
        eidtypeq,
        nxq, nyq, nzq, fsq % 面积分系数
    end
    
    methods
        function obj = mesh2d_fullquad(cell, varargin)
            obj = obj@ndg_lib.mesh.mesh2d(cell, varargin{:});
            % 计算每个面质量矩阵
            obj.Jq = obj.cell.proj_node2quad(obj.J);
            obj.invM = mass_matrix(obj);
            
            obj.rxq = obj.vol_scal(obj.rx); % 计算面积分转换系数
            obj.ryq = obj.vol_scal(obj.ry);
            obj.sxq = obj.vol_scal(obj.sx);
            obj.syq = obj.vol_scal(obj.sy);
            % 面积分 lift 矩阵
            [obj.nxq, obj.nyq, obj.nzq, obj.fsq] ...
                = face_scal(obj, obj.vx, obj.vy, obj.EToV);
            [obj.eidMq, obj.eidPq, obj.eidtypeq] = connect_quad_node(obj);
        end
    end
    
    methods(Access=private)
        function [nxq, nyq, nzq, fsq] = face_scal(obj, vx, vy, EToV)
            nxq = zeros(obj.cell.Nfq, obj.K);
            nyq = zeros(obj.cell.Nfq, obj.K);
            nzq = zeros(obj.cell.Nfq, obj.K);
            Nfq = obj.cell.Nfq/obj.cell.Nface; % 每个面 gauss 节点数
            % start index of each face node
            faceIndexStart = ones(obj.cell.Nface, 1); 
            for f = 2:obj.cell.Nface
                faceIndexStart(f) = faceIndexStart(f-1) + Nfq;
            end

            for f = 1:obj.cell.Nface
                Nfp = Nfq;
                face_x1 = vx(EToV(obj.cell.FToV(1,f), :))';
                face_x2 = vx(EToV(obj.cell.FToV(2,f), :))';
                face_y1 = vy(EToV(obj.cell.FToV(1,f), :))';
                face_y2 = vy(EToV(obj.cell.FToV(2,f), :))';

                ind = faceIndexStart(f):(faceIndexStart(f)+Nfp-1);
                nxq(ind, :) = repmat( (face_y2 - face_y1), Nfp, 1 );
                nyq(ind, :) = repmat(-(face_x2 - face_x1), Nfp, 1 );
            end

            % normalise
            Js = sqrt(nxq.*nxq+nyq.*nyq); 
            nxq = nxq./Js; 
            nyq = nyq./Js;
            fsq = bsxfun(@times, obj.cell.wbq, Js.*0.5);
        end
        
        function [rxq] = vol_scal(obj, rx)
            rxq = bsxfun(@times, obj.cell.wq, ...
                obj.cell.proj_node2quad(rx).*obj.Jq);
        end
    
        function invM = mass_matrix(obj)
            invM = zeros(obj.cell.Np, obj.cell.Np, obj.K);
            for k = 1:obj.K % 计算单元质量矩阵
                JVq = bsxfun(@times, obj.cell.wq.*obj.Jq(:, k), ...
                    obj.cell.Vq);
                mass_mat = obj.cell.Vq'*JVq;
                invM(:, :, k) = inv(mass_mat);
            end
        end
        
        function [eidMq, eidPq, eidtype] = connect_quad_node(obj)
            % 连接 gauss 积分节点
            
            % calculate the total nodes
            Nfnode = obj.cell.Nfq/obj.cell.Nface; % num nodes on each edge
            eidMq = zeros(obj.cell.Nfq, obj.K); 
            eidPq = zeros(obj.cell.Nfq, obj.K);
            eidtype = int8(zeros(obj.cell.Nfq, obj.K));
            
            sk = 1; list = 1:obj.cell.Nfq;
            for k1 = 1:obj.K
                xMq = obj.cell.proj_node2surf_quad(obj.x(:, k1));
                yMq = obj.cell.proj_node2surf_quad(obj.y(:, k1));
                for f1 = 1:obj.cell.Nface
                    k2 = obj.EToE(f1, k1);
                    xPq = obj.cell.proj_node2surf_quad(obj.x(:, k2));
                    yPq = obj.cell.proj_node2surf_quad(obj.y(:, k2));
                    type = int8(obj.EToBS(f1, k1));
                    for n = 1:Nfnode
                        loc_face_ind1 = (f1-1)*Nfnode + n;
                        eidMq(sk) = (k1-1)*obj.cell.Nfq + loc_face_ind1;

                        xpM = xMq(loc_face_ind1);
                        ypM = yMq(loc_face_ind1);

                        d12 = (xpM - xPq).^2 + (ypM - yPq).^2;
                        m = (d12 < 3e-16);
                        try
                            eidPq(sk) = (k2-1)*obj.cell.Nfq + list(m);
                        catch
                            keyboard;
                        end
                        eidtype(sk) = type;
                        sk = sk+1;
                    end
                end
            end% for
        end
    end
end

