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
        
        function obj = add_sponge(obj, vertlist)
            % 计算海绵层内单元距离边界长度
            obj. spg_delta = zeros(obj.cell.Np, obj.K);
            xb = obj.vx(vertlist);
            yb = obj.vy(vertlist);
            for k = 1:obj.K
                if obj.EToR(k) ~= ndg_lib.mesh_type.Sponge
                    continue;
                end
                
                for n = 1:obj.cell.Np
                    xi = obj.x(n, k);
                    yi = obj.y(n, k);
                    obj.spg_delta(n, k) = ...
                        min( sqrt( (xi - xb).^2 + (yi - yb).^2 ) );
                end
            end
        end
        
        function spg_sigma = cal_sponge_strength(obj, spg_len, max_sigma)
            % 计算海绵层内松弛系数 sigma
            spg_sigma = zeros(obj.cell.Np, obj.K);
            p = 3;
            if isempty(obj.spg_delta)
                return;
            end
            for k = 1:obj.K
                if obj.EToR(k) ~= ndg_lib.mesh_type.Sponge
                    continue;
                end
                
                spg_sigma(:, k) = ...
                    max_sigma*(1 - obj.spg_delta(:, k)/spg_len).^p;
            end
        end
        
        function draw(obj, varargin)
            switch nargin
                case 1
                    hold on;
                    patch('Vertices', [obj.vx(:), obj.vy(:)], ...
                        'Faces', obj.EToV', ...
                        'FaceColor', [0.8, 0.9, 1]);
                    plot(obj.x(:), obj.y(:), 'k.')
                case 2
                    hold on;
                    f = varargin{1};
                    EToV = ones(obj.K, 1)*obj.cell.Fmask(:)';
                    EToV = EToV + ( obj.cell.Np*(0:obj.K-1) )'...
                        *ones(1, obj.cell.Nfptotal);
                    patch(...
                        'Vertices', [obj.x(:), obj.y(:), f(:)], ...
                        'Faces', EToV, ...
                        'FaceColor', 'interp', ...
                        'FaceVertexCData', f(:));
                    plot(obj.x(:), obj.y(:), 'k.')
            end
        end% func
    end% methods
    
end

