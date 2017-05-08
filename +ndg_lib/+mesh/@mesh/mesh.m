classdef mesh
    %@STD_MESH Summary of this class goes here
    %   Detailed explanation goes here
    
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
    properties(Abstract, SetAccess=private)
        EToE, EToF
        x, y, z
        rx, ry, rz
        sx, sy, sz
        tx, ty, tz
        J
    end
    
    % elemental edge information
    properties(Abstract, SetAccess=private)
        nx, ny, nz
        Js
    end
    
    % edge information
    properties(Abstract, SetAccess=private)
        Nedge
        kM, kP
        fM, fP
        ftype
        idM, idP
        fpM, fpP
        fscal
        fnxM, fnyM, fnzM
    end
    
    % I/O methods
    methods(Abstract)
        % read mesh information from input file, I/O methods
        [Nv, vx, vy, vz, K, EToV, EToR, EToBS] = read_from_file(casename)
    end
    
    % Abstract methods
    methods(Abstract, Access = protected)
        ele_connect(obj) % get EToE and EToF through EToV
        ele_vol_factor(obj) % get volume infomation () of each element
        ele_suf_factor(obj) % get out normal vector of each elemental edges
        edge_connect(obj) % get edge connection through element connection   
    end
    
    % public methods
    methods
        function nodeQ = proj_vert2node(obj, vertQ)
            % project scalars from mesh verts to nodes
            ele_vQ = vertQ(obj.EToV); 
            nodeQ = obj.cell.project_vert2node(ele_vQ);
        end
    end
    
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
        end% func
        
        function [x, y, z] = ele_node_project(obj, vx, vy, vz)
            x = obj.proj_vert2node(vx);
            y = obj.proj_vert2node(vy);
            z = obj.proj_vert2node(vz);
        end% func
        
    end
    
end

