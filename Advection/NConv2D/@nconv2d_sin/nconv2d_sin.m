classdef nconv2d_sin < nconv2d
    %NCONV2D_SIN Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = nconv2d_sin(varargin)
            if( isa(varargin{2}, 'char') )
                N = varargin{1};
                casename = varargin{2};
                cell_type = varargin{3};
                input_type = 'file';
                input_var = {N, cell_type, casename};
            elseif( isa(varargin{2}, 'double') )
                N = varargin{1}; % order of the basis function
                M = varargin{2}; % # of elements
                cell_type = varargin{3}; % cell type
                xlim = [0, 2*pi]; ylim = [0, 2*pi]; % computetion domain
                Mx = M; My = M;
                % set conditions for the sourth, north, west and east
                % boundaries
                cl_bc = ndg_lib.bc_type.Clamped;
                bc_type = [cl_bc, cl_bc, cl_bc, cl_bc]; 

                input_type = 'uniform';
                input_var = {N, cell_type, xlim, ylim, Mx, My, bc_type};
            end
            obj = obj@nconv2d(input_type, input_var);
            obj.ftime = 0.6;
            obj.init();
        end% func
        
        function init(obj)
            obj.f_Q = obj.init_func( obj.mesh.x, obj.mesh.y );
            obj.f_extQ = obj.f_Q;
        end% func
    end
    
    methods(Access=protected)   
        
        function [ f_extQ ] = ext_func(obj, time)
            % calculate the exact solution based on the character method
            xn = obj.mesh.x; % last step position
            yn = obj.mesh.y;
            eof = 1;
            while ( any(eof > 1e-10) )
                u0 = obj.init_func(xn, yn);
                x0 = obj.mesh.x - u0*time; % postion in next step
                y0 = obj.mesh.y - u0*time;
                eof = sqrt( (x0(:) - xn(:)).^2 + (y0(:) - yn(:)).^2 );
                xn = x0; 
                yn = y0;
            end
            
            f_extQ = obj.init_func(xn, yn);
        end% func
    end
    
    methods(Static)
        function [ f_init ] = init_func(x, y)
            f_init = (1 - cos( x )).*(1 - cos( y ))/2;
        end% func
    end
    
end

