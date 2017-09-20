%> @brief Create a generalized 2D unstructured mesh class for triangle or
%> quadrilateral meshes.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDG-FEM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef UMeshUnion2d < UMeshUnion
    
    properties(Hidden=true)
        point_h
    end
    
    methods(Hidden, Access=protected)
        [rx, ry, rz, sx, sy, sz, tx, ty, tz, J] = getElementalNodeInfo(obj)
        
        function [ uedge ] = setUEdgeClass(obj, BCToV)
            uedge = UEdgeLine(obj, BCToV);
        end
    end% methods
    
    %% public methods
    methods(Access=public)
        function obj = UMeshUnion2d(cell, varargin)
            
            switch varargin{1}
                case 'file'
                    [Nv, vx, vy, K, EToV, EToR, BCToV] ...
                        = read_from_file( varargin{2} );
                case 'uniform'
                    [Nv, vx, vy, K, EToV, EToR, BCToV] ...
                        = create_unfirm_mesh(cell, varargin{2});
                case 'variable'
                    [Nv, vx, vy, K, EToV, EToR, BCToV] ...
                        = check_var_input(cell, varargin{2});
                otherwise
            end
               
            vz = zeros(size(vx)); % vz is all zeros
            obj = obj@UMeshUnion(cell, Nv, vx, vy, vz, K, EToV, EToR, BCToV);
        end% func
        
%         function obj = add_sponge(obj, vertlist)
%             % Calculate the distance from the boundary for each nodes
%             obj. spg_delta = zeros(obj.cell.Np, obj.K);
%             xb = obj.vx(vertlist);
%             yb = obj.vy(vertlist);
%             for k = 1:obj.K
%                 if obj.EToR(k) ~= ndg_lib.mesh_type.Sponge
%                     continue;
%                 end
%                 
%                 for n = 1:obj.cell.Np
%                     xi = obj.x(n, k);
%                     yi = obj.y(n, k);
%                     obj.spg_delta(n, k) = ...
%                         min( sqrt( (xi - xb).^2 + (yi - yb).^2 ) );
%                 end
%             end
%         end
        
%         function spg_sigma = cal_sponge_strength(obj, spg_len, max_sigma)
%             % 计算海绵层内松弛系数 sigma
%             spg_sigma = zeros(obj.cell.Np, obj.K);
%             p = 3;
%             if isempty(obj.spg_delta)
%                 return;
%             end
%             for k = 1:obj.K
%                 if obj.EToR(k) ~= ndg_lib.mesh_type.Sponge
%                     continue;
%                 end
%                 
%                 spg_sigma(:, k) = ...
%                     max_sigma*(1 - obj.spg_delta(:, k)/spg_len).^p;
%             end
%         end
        
        function draw(obj, varargin)
            switch nargin
                case 1
                    hold on;
                    if ( isempty(obj.draw_h) || ~isvalid(obj.draw_h))
                        obj.draw_h = patch(...
                            'Vertices', [obj.vx(:), obj.vy(:)], ...
                            'Faces', obj.EToV', ...
                            'FaceColor', [0.8, 0.9, 1]);
                        obj.point_h = plot(obj.x(:), obj.y(:), 'k.');
                        box on;
                    else
                        set(obj.draw_h, 'Vertices', [obj.x(:), obj.y(:)]);
                    end% if
                case 2
                    hold on;
                    f = varargin{1};
                    if ( isempty(obj.draw_h) || ~isvalid(obj.draw_h))
                        EToV = ones(obj.K, 1)*obj.cell.Fmask(:)';
                        EToV = EToV + ( obj.cell.Np*(0:obj.K-1) )'...
                            *ones(1, obj.cell.Nfptotal);
                        obj.draw_h = patch(...
                            'Vertices', [obj.x(:), obj.y(:), f(:)], ...
                            'Faces', EToV, ...
                            'FaceColor', 'interp', ...
                            'FaceVertexCData', f(:));
                        obj.point_h = plot(obj.x(:), obj.y(:), 'k.');
                        box on;
                    else % if the figure exists
                        set(obj.draw_h, ...
                            'Vertices', [obj.x(:), obj.y(:), f(:)],...
                            'FaceVertexCData', f(:));
                    end
            end
        end% func
    end% methods
    
end

function [Nv, vx, vy, K, EToV, EToR, BCToV] = check_var_input(cell, input)
% check the input variables for initlizing the mesh object.
Nv = input{1};
vx = input{2}; 
vy = input{3}; 
K  = input{4}; 
EToV = input{5};
EToR = int8(input{6}); 
BCToV = int8(input{7});

% check the number of the vertex
if ( numel(vx) ~= Nv ) || ( numel(vy) ~= Nv )
    error(['The length of input vertex coordinate "vx" or "vy" ', ...
        'is not equal to Nv']);
end% if

% check the number of the elements
if ( size(EToV, 2) ~= K ) ||( numel(EToR) ~= K )
    error(['The length of input "EToV" or "EToR" ', ...
        'is not equal to K']);
end% if

if ( size(EToV, 1) ~= cell.Nv) || ( size(BCToV, 1) ~= 3 )
    error('The numbers of vertex in "EToV" and "EToBS" is not fault');
end% if
end% func

function [Nv, vx, vy, K, EToV, EToR, EToBS] = create_unfirm_mesh(cell, var)
% Generate the informations in the unifrom mesh.
xlim = var{1}; 
ylim = var{2};
Mx = var{3}; 
My = var{4};
facetype = var{5};

switch cell.type
    case StdCellType.Tri % uniform triangle mesh
        [Nv, vx, vy, K, EToV, EToR, EToBS] ...
            = tri_uniform_mesh(xlim, ylim, Mx, My, facetype);
    case StdCellType.Quad % uniform quadrilateral mesh
        [Nv, vx, vy, K, EToV, EToR, EToBS] ...
            = quad_uniform_mesh(xlim, ylim, Mx, My, facetype);
end% switch

end% func
