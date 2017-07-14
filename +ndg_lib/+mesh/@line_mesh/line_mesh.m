classdef line_mesh < ndg_lib.mesh.mesh
    %MESH_LINE Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(Static)
        [Nv, vx, K, EToV, EToR, EToBS] = read_from_file(casename)
    end
    
    methods(Hidden, Access = protected)
        [rx, ry, rz, sx, sy, sz, tx, ty, tz, J] = ele_vol_factor(obj)
        [nx, ny, nz, Js] = ele_suf_factor(obj, vx, vy, vz, EToV)
        
        function Eind = get_Eind(obj)
            Eind = zeros(obj.cell.Nface, obj.K);
            for f = 1:obj.cell.Nface
                Eind(f, :) = obj.EToV(obj.cell.FToV(1,f), :);
            end
        end% func
    end% methods
    
    methods
        obj = refine(obj, refine_level);
        
        function obj = line_mesh(cell, varargin)
            if(cell.type ~= ndg_lib.std_cell_type.Line) % check input cell
                error(['Input cell type ', cell.type, 'is not line!'])
            end
            if(nargin == 2) % case input
                casename = varargin{1};
                [Nv, vx, K, EToV, EToR, EToBS] = read_from_file(casename);
            elseif(nargin == 7) % parameters input
                Nv = varargin{1};
                vx = varargin{2};
                K = varargin{3};
                EToV = varargin{4};
                EToR = varargin{5};
                EToBS = varargin{6};
            else
                error(['The number of input parameters ', num2str(nargin), ...
                    ' is incorrect!']);
            end
            vy = zeros(size(vx)); % vy is all zeros
            vz = zeros(size(vx)); % vz is all zeros
            
            obj = obj@ndg_lib.mesh.mesh(cell, ...
                Nv, vx, vy, vz, K, EToV, EToR, EToBS);
        end% func
        
        function draw(obj)
            plot(obj.x, zeros(obj.cell.Np, obj.K), '.-');
        end
    end
    
end

