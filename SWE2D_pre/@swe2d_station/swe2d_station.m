classdef swe2d_station < swe2d
    %SWE2D_STATION The station problem to test the well-balanced problem.
    %   Detailed explanation of this problem could be refered to
    %   Xing et al. (Advances in Water Resources, 2010).
    
    properties(Constant)
        hmin = 1e-4;
    end
    
    methods
        function obj = swe2d_station(varargin)
            [ input ] = set_case_parameter( varargin{:} );
            obj = obj@swe2d(input{:});
            obj.ftime = 10;
            obj.init();
        end% func
        
        function init( obj )
            r2 = ( (obj.mesh.x - 0.5).^2 + (obj.mesh.y - 0.5).^2 );
            obj.bot = max(0, 0.25 - 5*r2);
            %obj.bot = 5*r2 - 0.5;
            h = 0.2 - obj.bot;
            h( h < 0 ) = 0;
            
            obj.f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            obj.f_Q(:, :, 1) = h;
            obj.f_Q(:, :, 2:3) = 0;
        end% func
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
    M = varargin{2}; % # of elements on each axis
    type = varargin{3};
    xlim = [0, 1]; 
    ylim = [0, 1];
    Mx = M; 
    My = M;
    zg_bc = ndg_lib.bc_type.ZeroGrad;
    bc_type = [zg_bc, zg_bc, zg_bc, zg_bc];
    input = {N, xlim, ylim, Mx, My, type, bc_type};
end% switch

end% func

