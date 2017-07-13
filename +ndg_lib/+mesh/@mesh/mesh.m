classdef mesh
    %@STD_MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    %%
    properties
        cell
        K
        Nv
        EToV % element to vertex information
        EToR@int8
        EToBS
        vx, vy, vz
        J
    end
    
    % elemental volume infomation
    properties(SetAccess=protected)
        EToE, EToF
        x, y, z
        rx, ry, rz
        sx, sy, sz
        tx, ty, tz
        vol
        spg_delta % 海绵层距离边界点距离
    end
    
    properties(SetAccess=protected)
        eidM, eidP
        eidtype@int8
        eidfscal
    end
    
    % elemental edge information
    properties(SetAccess=protected)
        nx, ny, nz
        Js, elen
    end
    
    % edge information
    properties(SetAccess=protected)
        Nedge
        Nnode
        kM, kP
        fM, fP
        ftype@int8
        idM, idP
        fpM, fpP
        fscal
        fnxM, fnyM, fnzM
    end
    
    
    %% Abstract methods
    methods(Abstract, Access = protected)
        ele_connect(obj) % get EToE and EToF through EToV
        ele_vol_factor(obj) % get volume infomation () of each element
        ele_suf_factor(obj) % get out normal vector of each elemental edges
        edge_connect(obj) % get edge connection through element connection   
    end
    methods(Abstract)
        draw(obj); 
    end
    
    %% public methods
    methods
        function nodeQ = proj_vert2node(obj, vertQ)
            % project scalars from mesh verts to nodes
            ele_vQ = vertQ(obj.EToV); 
            nodeQ = obj.cell.project_vert2node(ele_vQ);
        end
        
        function c_m = cell_mean(obj, f_Q)
            % calculate the cell averaged values
%             par = bsxfun(@times, obj.cell.w, obj.J);
%             c_m = sum( par.*f_Q )./sum(par);
            [ c_m ] = cell_mean( f_Q, obj.cell.w, obj.J );
        end
        
        function f_m = face_mean(obj, f_Q)
            % calculate the averaged values of each faces
            par = bsxfun(@times, obj.cell.ws, obj.Js);
            f_m = zeros(obj.cell.Nface, obj.K);
            f_M = f_Q(obj.eidM);
            tmp = par.*f_M;
            st = 1;
            if (obj.cell.Nfp == 1)
                f_m = tmp;
            else
                for f = 1:obj.cell.Nface
                    sk = st + obj.cell.Nfp(f);
                    row = st:(sk-1);
                    f_m(f, :) = sum( tmp(row, :) )./sum(par(row, :));
                    st = sk;
                end
            end
        end
    end
    
    % I/O methods
    methods(Abstract)
        % read mesh information from input file, I/O methods
        [Nv, vx, vy, vz, K, EToV, EToR, EToBS] = read_from_file(casename)
    end
    
    %% 公有函数
    methods
        function obj = mesh(cell, Nv, vx, vy, vz, K, EToV, EToR, EToBS)
            obj.cell = cell; obj.Nv = Nv;
            obj.vx = vx; obj.vy = vy; obj.vz = vz;
            obj.K  = K; obj.EToV = EToV;
            obj.EToR = int8(EToR); obj.EToBS = int8(EToBS);
            
            [obj.EToE, obj.EToF] = ele_connect(obj, obj.EToV);
            [obj.x, obj.y, obj.z] ...
                = ele_node_project(obj, obj.vx, obj.vy, obj.vz);
            
            [obj.rx, obj.ry, obj.rz, ...
                obj.sx, obj.sy, obj.sz, ...
                obj.tx, obj.ty, obj.tz, obj.J] = ele_vol_factor(obj);
            
            [obj.nx, obj.ny, obj.nz, obj.Js] ...
                = ele_suf_factor(obj, obj.vx, obj.vy, obj.EToV);
            
            [obj.Nedge, obj.Nnode, ...
                obj.kM, obj.kP, obj.fM, obj.fP, obj.ftype, ...
                obj.idM, obj.idP, obj.fpM, obj.fpP, obj.fscal, ...
                obj.fnxM, obj.fnyM, obj.fnzM] ...
                = edge_connect(obj, obj.EToV, obj.EToE, obj.EToF, obj.EToBS);
            
            [obj.eidM, obj.eidP, obj.eidtype] = ele_node_connect(obj);
            
            [obj.vol, obj.elen] = ele_scale(obj);
        end% func
    end
    
    %% 私有函数
    methods(Access = protected)
        
        function [vol, edge_len] = ele_scale(obj)
            vol = sum( bsxfun(@times, obj.cell.w, obj.J) );
            
            par = bsxfun(@times, obj.cell.ws, obj.Js);
            edge_len = zeros(obj.cell.Nface, obj.K);
            if (obj.cell.Nfp == 1)
                edge_len = par;
            else
                st = 1;
                for f = 1:obj.cell.Nface
                    sk = st + obj.cell.Nfp(f);
                    row = st:(sk-1);
                    edge_len(f, :) = sum(par(row, :));
                    st = sk;
                end
            end
        end% func
        
        function [x, y, z] = ele_node_project(obj, vx, vy, vz)
            % 通过单元顶点坐标获取单元内节点坐标
            x = obj.proj_vert2node(vx);
            y = obj.proj_vert2node(vy);
            z = obj.proj_vert2node(vz);
        end% func
        
        function [eidM, eidP, eidtype] = ele_node_connect(obj)
            % 连接相邻边界节点编号
            Nfptotal = obj.cell.Nfptotal;
            eidM = zeros(Nfptotal, obj.K);
            eidP = zeros(Nfptotal, obj.K);
            eidtype = int8(zeros(Nfptotal, obj.K));
            
            Np = obj.cell.Np;
            for k1 = 1:obj.K
                
                eidM(:, k1) = (k1-1)*Np + obj.cell.Fmask(:);
                f1 = 1;
                nfp1 = 1;
                for n = 1:Nfptotal
                    if nfp1 > obj.cell.Nfp(f1) % 若面节点编号超过f1面节点总数
                        nfp1 = 1;
                        f1 = f1 + 1; % 面编号+1
                    end
                    eidtype(n, k1) = int8(obj.EToBS(f1, k1));
                    
                    k2 = obj.EToE(f1, k1); % 相邻单元
                    f2 = obj.EToF(f1, k1); % 相邻面
                    
                    xpM = obj.x(eidM(n, k1)); % 本单元边界节点坐标
                    ypM = obj.y(eidM(n, k1));
                    zpM = obj.z(eidM(n, k1));
                    % 相邻面上所有节点坐标
                    ind2 = (k2-1)*Np + obj.cell.Fmask(:, f2); 
                    xP = obj.x( ind2 );
                    yP = obj.y( ind2 );
                    zP = obj.z( ind2 );
                    d12 = (xpM - xP).^2 + (ypM - yP).^2 + (zpM - zP).^2;
                    m = (d12 < 3e-16);
%                     try
                    eidP(n, k1) = (k2-1)*Np + obj.cell.Fmask(m, f2);
%                     catch
%                         keyboard
%                     end
                    nfp1 = nfp1 + 1;
                end
            end
        end
    end
    
end

