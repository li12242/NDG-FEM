classdef mesh2d < ndg_lib.mesh.mesh
    %MESH2D 2d mesh object
    %   Create the 2d mesh object. The object is inherited from the mesh
    %   object and has special treatment for the 2d elements (triangle and
    %   quadrilateral). The mesh object contains only one type of elements,
    %   and allocate multiple choices to create the object. The input
    %   methods include
    %   * use file;
    %   * use parameters to create uniform mesh;
    %   * use mesh variables to create the mesh;
    %

    methods(Static)
        [Nv, vx, vy, K, EToV, EToR, EToBS] = read_from_file(casename)
    end
    methods(Abstract, Static)
        [Nv, vx, vy, K, EToV, EToR, EToBS] ...
            = uniform_mesh(xlim, ylim, Mx, My, facetype)
    end
    
    %% private methods
    methods(Hidden, Access=protected)
        [rx, ry, rz, sx, sy, sz, tx, ty, tz, J] = ele_vol_factor(obj)
        [nx, ny, nz, Js] = ele_suf_factor(obj, vx, vy, vz, EToV)
        
        function Eind = get_Eind(obj)
            Eind = zeros(obj.cell.Nface, obj.K);
            for f = 1:obj.cell.Nface
                v1 = obj.EToV(obj.cell.FToV(1,f), :);
                v2 = obj.EToV(obj.cell.FToV(2,f), :);
                % calculate the indicator for each edge
                Eind(f, :) = min(v1, v2)*obj.Nv + max(v1, v2);
            end
        end
    end% methods
    
    %% public methods
    methods(Access=public)
        function obj = mesh2d(cell, varargin)
            
            switch varargin{1}
                case 'file'
                    [Nv, vx, vy, K, EToV, EToR, EToBS] ...
                        = read_from_file( varargin{2} );
                case 'uniform'
                    [Nv, vx, vy, K, EToV, EToR, EToBS] ...
                        = unfirm_mesh(cell, varargin{2});
                case 'variable'
                    [Nv, vx, vy, K, EToV, EToR, EToBS] ...
                        = check_var_input(cell, varargin{2});
                otherwise
            end
               
            vz = zeros(size(vx)); % vz is all zeros
            obj = obj@ndg_lib.mesh.mesh(cell, ...
                Nv, vx, vy, vz, K, EToV, EToR, EToBS);
        end% func
        
        function obj = add_sponge(obj, vertlist)
            % Calculate the distance from the boundary for each nodes
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

function [Nv, vx, vy, K, EToV, EToR, EToBS] = check_var_input(cell, input)
% check the input variables for initlizing the mesh object.
Nv = input{1};
vx = input{2}; 
vy = input{3}; 
K  = input{4}; 
EToV = input{5};
EToR = int8(input{6}); 
EToBS = int8(input{7});

% check the # of the vertex
if ( numel(vx) ~= Nv ) || ( numel(vy) ~= Nv )
    error(['The length of input vertex coordinate "vx" or "vy" ', ...
        'is not equal to Nv']);
end% func

% check the # of the elements
if ( size(EToV, 2) ~= K )||( size(EToBS, 2) ~= K ) ||( numel(EToR) ~= K )
    error(['The length of input "EToV", "EToR" or "EToBS" ', ...
        'is not equal to K']);
end% func

if ( size(EToV, 1) ~= cell.Nv) || ( size(EToBS, 1) ~= cell.Nv )
    error('The numbers of vertex in "EToV" and "EToBS" is not fault');
end% func
end% func

function [Nv, vx, vy, K, EToV, EToR, EToBS] = unfirm_mesh(cell, var)
% Generate the informations in the unifrom mesh.
xlim = var{1}; 
ylim = var{2};
Mx = var{3}; 
My = var{4};
facetype = var{5};

switch cell.type
    case ndg_lib.std_cell_type.Tri % uniform triangle mesh
        [Nv, vx, vy, K, EToV, EToR, EToBS] ...
            = tri_uniform_mesh(xlim, ylim, Mx, My, facetype);
    case ndg_lib.std_cell_type.Quad % uniform quadrilateral mesh
        [Nv, vx, vy, K, EToV, EToR, EToBS] ...
            = quad_uniform_mesh(xlim, ylim, Mx, My, facetype);
end% switch

end% func
