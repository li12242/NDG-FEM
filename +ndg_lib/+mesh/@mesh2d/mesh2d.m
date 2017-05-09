classdef mesh2d < ndg_lib.mesh.mesh
    %@STD_MESH2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cell
        K
        Nv
        EToV
        EToR
        EToBS
        vx, vy, vz
    end
    
    % elemental infomation
    properties(SetAccess=private)
        EToE, EToF
        x, y, z
        rx, ry, rz
        sx, sy, sz
        tx, ty, tz
        nx, ny, nz
        J
        Js
    end
    
    properties(SetAccess=private)
        eidM, eidP
        eidtype@uint8
        eidfscal
    end
    
    properties(SetAccess=private)
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
    
    methods(Static)
        [Nv, vx, vy, K, EToV, EToR, EToBS] = read_from_file(casename)
    end
    methods(Access = protected)
        [EToE, EToF] = ele_connect(obj, EToV)
        [rx, ry, rz, sx, sy, sz, tx, ty, tz, J] = ele_vol_factor(obj)
        [nx, ny, nz, Js] = ele_suf_factor(obj, vx, vy, EToV)
        [Nedge, Nnode, kM, kP, fM, fP, ftype, ...
            idM, idP, fpM, fpP, fscal, fnxM, fnyM, fnzM] = ...
            edge_connect(obj, EToV, EToE, EToF, EToBS)
    end% methods
    
    methods
        function obj = mesh2d(cell, varargin)
            
            if(nargin == 2) % case input
                casename = varargin{1};
                [Nv, vx, vy, K, EToV, EToR, EToBS] = read_from_file(casename);
            elseif(nargin == 8) % parameters input
                Nv = varargin{1};
                vx = varargin{2};
                vy = varargin{3};
                K = varargin{4};
                EToV = varargin{5};
                EToR = varargin{6};
                EToBS = varargin{7};
            else
                error(['The number of input parameters ', num2str(nargin), ...
                    ' is incorrect!']);
            end
            vz = zeros(size(vx)); % vz is all zeros
            
            obj = obj@ndg_lib.mesh.mesh(cell, ...
                Nv, vx, vy, vz, K, EToV, EToR, EToBS);
                        
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
        end% func
    end% methods
    
end

