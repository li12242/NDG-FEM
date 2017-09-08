classdef swe2d_dbw < swe2d
    %SWE2D_DBW The dam break problem with wet condition at the downstream.
    %   
    
    properties(Constant)
        dam_pos = 500
        h0 = 10 % water depth at the upstream (left)
        h1 = 2 % water depth at the downstream (right)
    end% properties
    
    properties
        theta = pi/6; % clockwise rotation angle
    end% properties
    
    properties(Constant)
        hmin = 1e-4 % the threshold of the water depth
    end% properties
    
    methods(Access=protected)
        [ f_ext ] = ext_func( obj, t )
    end% methods
    
    methods
        function [ f_ext ] = get_ext(obj, t)
            [ f_ext ] = ext_func( obj, t );
        end
        
        function obj = swe2d_dbw(varargin)
            [ input ] = set_case_parameter( varargin{:} );
            obj = obj@swe2d(input{:});
            obj.ftime = 20;
            obj.init();
        end% func
        
        function init(obj)
            h = obj.h1*ones(obj.mesh.cell.Np, obj.mesh.K);
            % the inverse rotation matrix
            T = [cos(-obj.theta), sin(-obj.theta); 
                -sin(-obj.theta), cos(-obj.theta)]; 

            x = zeros(size(obj.mesh.x));
            temp = T*[ obj.mesh.x(:)'; obj.mesh.y(:)'; ];
            x(:) = temp(1, :);
            xc = obj.mesh.cell_mean(x);
            h(:, xc<obj.dam_pos) = obj.h0;
            
            obj.f_Q(:, :, 1) = h;
            obj.f_Q(:, :, 2:3) = 0;
            obj.bot = zeros(obj.mesh.cell.Np, obj.mesh.K); % bottom
        end% func
        
        result(obj);
        [ obj ] = SL_RK45( obj );
    end
    
end

function [ input ] = set_case_parameter( varargin )
% set the range of the computation domain and the boundary conditions
if( isa(varargin{2}, 'char') )
    N = varargin{1};
    casename = varargin{2};
    type = varargin{3};
    input = {N, casename, type};
elseif( isa(varargin{2}, 'double') )
    N = varargin{1};
    M = varargin{2};
    type = varargin{3};
    xlim = [0, 1000]; 
    ylim = [0, 20];
    Mx = M; My = 1;
    zg_bc = ndg_lib.bc_type.ZeroGrad;
    bc_type = [zg_bc, zg_bc, zg_bc, zg_bc];
    input = {N, xlim, ylim, Mx, My, type, bc_type};
end% switch

end% func

