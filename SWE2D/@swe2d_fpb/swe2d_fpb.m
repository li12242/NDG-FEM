classdef swe2d_fpb < swe2d
    %SWE2D_FPB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        hmin = 1e-4
        h0 = 10
        a = 3e3
        B = 5
        k = 0.002 % 线性摩阻系数
        T = 2*672
    end
    
    methods(Access=protected)
        [ f_ext ] = extval(obj, time);
    end
    
    methods
        function obj = swe2d_fpb(varargin)
            switch nargin
                case 1
                    mesh = varargin(1);
                    if ( ~isa(mesh, 'ndg_lib.mesh.tri_mesh') || ...
                            ~isa(mesh, 'ndg_lib.mesh.quad_mesh') )
                        error(['The input is not a triangle or ',... 
                            'quadrilateral mesh object!']);
                    end
                case 3
                    N = varargin{1}; 
                    M = varargin{2}; 
                    type = varargin{3};
                    mesh = uniform_mesh( N, M, type );
                otherwise
                    error('The number of input variable is incorrect.');
            end% switch
            
            obj = obj@swe2d(mesh);
            obj.cfl = 0.2;
            obj.ftime = 4*obj.T;
            obj.init();
        end
        
        function init(obj)
            r2 = obj.mesh.x.^2 + obj.mesh.y.^2;
            obj.bot = r2*(obj.h0/obj.a^2);
            obj.f_Q = obj.extval(0);
        end
    end
    
end

