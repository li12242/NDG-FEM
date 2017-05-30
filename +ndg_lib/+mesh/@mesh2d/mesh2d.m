classdef mesh2d < ndg_lib.mesh.mesh
    %@STD_MESH2D Summary of this class goes here
    %   Detailed explanation goes here

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
        end% func
        
        function draw(obj)
            hold on;
            patch('Vertices', [obj.vx(:), obj.vy(:)], ...
                'Faces', obj.EToV', ...
                'FaceColor', [0.8, 0.9, 1]);
            plot(obj.x(:), obj.y(:), 'k.')
        end% func
    end% methods
    
end

