classdef mesh
    %@STD_MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    %%
    properties(Abstract)
        cell
        K
        Nv
        EToV % element to vertex information
        EToR
        EToBS
        vx, vy, vz
    end
    
    % elemental volume infomation
    properties(SetAccess=protected)
        EToE, EToF
        x, y, z
        rx, ry, rz
        sx, sy, sz
        tx, ty, tz
        J, vol
    end
    
    properties(SetAccess=protected)
        eidM, eidP
        eidtype@uint8
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
        ftype@uint8
        idM, idP
        fpM, fpP
        fscal
        fnxM, fnyM, fnzM
    end
    
    %% I/O methods
    methods(Abstract)
        % read mesh information from input file, I/O methods
        [Nv, vx, vy, vz, K, EToV, EToR, EToBS] = read_from_file(casename)
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
            par = bsxfun(@times, obj.cell.w, obj.J);
            c_m = sum( par.*f_Q )./sum(par);
        end
        
        function f_m = face_mean(obj, f_Q)
            % calculate the averaged values of each faces
            par = bsxfun(@times, obj.cell.ws, obj.Js);
            f_m = zeros(obj.cell.Nface, 1);
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
    
    %% private methods
    methods(Access = protected)
        function obj = mesh(cell, Nv, vx, vy, vz, K, EToV, EToR, EToBS)
            obj.cell = cell;
            obj.Nv = Nv;
            obj.vx = vx;
            obj.vy = vy;
            obj.vz = vz;
            obj.K  = K;
            obj.EToV = EToV;
            obj.EToR = EToR;
            obj.EToBS = EToBS;
            
            [obj.EToE, obj.EToF] = ele_connect(obj, obj.EToV);
            [obj.x, obj.y, obj.z] = ele_node_project(obj, obj.vx, obj.vy, obj.vz);
            [obj.rx, obj.ry, obj.rz, obj.sx, obj.sy, obj.sz, ...
                obj.tx, obj.ty, obj.tz, obj.J] = ele_vol_factor(obj);
            [obj.nx, obj.ny, obj.nz, obj.Js] = ...
                ele_suf_factor(obj, obj.vx, obj.vy, obj.EToV);
            [obj.Nedge, obj.Nnode, obj.kM, obj.kP, obj.fM, obj.fP, obj.ftype, ...
                obj.idM, obj.idP, obj.fpM, obj.fpP, obj.fscal, ...
                obj.fnxM, obj.fnyM, obj.fnzM] = ...
                edge_connect(obj, obj.EToV, obj.EToE, obj.EToF, obj.EToBS);
            [obj.eidM, obj.eidP, obj.eidtype, obj.eidfscal] = ele_suf_connect(obj);
            [obj.vol, obj.elen] = ele_scale(obj);
        end% func
        
        function [vol, edge_len] = ele_scale(obj)
            vol = sum( bsxfun(@times, obj.cell.w, obj.J) );
            
            par = bsxfun(@times, obj.cell.ws, obj.Js);
            edge_len = zeros(obj.cell.Nface, 1);
            if (obj.cell.Nfp == 1)
                edge_len = par;
            else
                for f = 1:obj.cell.Nface
                    sk = st + obj.cell.Nfp(f);
                    row = st:(sk-1);
                    edge_len(f, :) = sum(par(row, :));
                    st = sk;
                end
            end
        end% func
        
        function [x, y, z] = ele_node_project(obj, vx, vy, vz)
            x = obj.proj_vert2node(vx);
            y = obj.proj_vert2node(vy);
            z = obj.proj_vert2node(vz);
        end% func
        
        function [eidM, eidP, eidtype, eidfscal] = ele_suf_connect(obj)
            Nfptotal = obj.cell.Nfptotal;
            eidM = zeros(Nfptotal, obj.K);
            eidP = zeros(Nfptotal, obj.K);
            eidtype = uint8(zeros(Nfptotal, obj.K));
            eidfscal = zeros(Nfptotal, obj.K);
            faceIndexStart = ones(obj.cell.Nface, 1); % start index of each face node
            for f = 2:obj.cell.Nface
                faceIndexStart(f) = faceIndexStart(f-1) + obj.cell.Nfp(f-1);
            end
            sk = 1;
            for f = 1:obj.Nedge
                Nfp = obj.cell.Nfp(obj.fM(f));
                k1 = obj.kM(f);
                list = sk:(sk+Nfp-1);
                eidM(obj.fpM(list), k1) = obj.idM(list);
                eidP(obj.fpM(list), k1) = obj.idP(list);
                eidtype(obj.fpM(list), k1) = obj.ftype(f);
                eidfscal(obj.fpM(list), k1) = obj.fscal(list);
                % adjacent element
                k2 = obj.kP(f);
                eidM(obj.fpP(list), k2) = obj.idP(list);
                eidP(obj.fpP(list), k2) = obj.idM(list);
                eidtype(obj.fpP(list), k2) = obj.ftype(f);
                eidfscal(obj.fpP(list), k2) = obj.fscal(list);
                sk = sk + Nfp;
            end
        end
    end
    
end

