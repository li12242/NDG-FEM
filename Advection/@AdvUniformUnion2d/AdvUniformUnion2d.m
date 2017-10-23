classdef AdvUniformUnion2d < handle
    %ADVUNIFORMUNION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        %> initial position
        x0 = -0.5
        %> initial position
        y0 = -0.5
        %> const flow field
        u0 = 0.5
        %> const flow field
        v0 = 0.5
    end
    
    properties(SetAccess = protected )
        %> physical field
        phys
        %> flow velocity
        u
        %> flow velocity
        v
    end
    
    properties
        %> scal
        fvar
    end
    
    methods
        function solve(obj)
            obj.fvar = AdvSolver2d(obj);
        end
        
        function obj = AdvUniformUnion2d(N, M, type)
            % check input and create uniform mesh
            mesh = makeUniformMesh(N, M, type);            
            obj.phys = NdgPhys(1, mesh);
            % set initial condition
            obj.init;
        end% func
        
        function init(obj)
            obj.phys.finalTime = 2;
            mesh = obj.phys.mesh;
            vel = sqrt( obj.u0.^2 + obj.v0^2 );
            range = max( mesh.x(:, 1) ) - min( mesh.x(:, 1) );
            dt = range/vel/(2*mesh.cell.N + 1);
            
            obj.phys.setTimeIntervalType( NdgIntervalType.Constant,  dt);
            obj.phys.setLimiterType( NdgLimiterType.None );
            obj.phys.setTemporalDiscreteType( NdgTemporalDiscreteType.RK45 );
            obj.phys.setBoundaryType( NdgBCType.None );
            obj.phys.setOutputFile( 'advUniformMesh', ...
                NdgIntervalType.DeltaTime, 0.002 );
            
            K = mesh.K;
            Np = mesh.cell.Np;
            obj.u = obj.u0 * ones(Np, K);
            obj.v = obj.v0 * ones(Np, K);
            obj.fvar = obj.ext_func(0);
        end% func
    end
    
    methods(Access = protected)
        %> the exact function
        function f_ext = ext_func(obj, time)
            xc = obj.x0 + obj.u.*time;
            yc = obj.y0 + obj.v.*time;
            
            mesh = obj.phys.mesh;
            sigma = 125*1e3/(33*33);
            t = -( (mesh.x-xc).^2+(mesh.y-yc).^2 )*sigma;
            f_ext = exp(t);
        end% func
    end
    
end

function mesh = makeUniformMesh(N, M, type)
bctype = [NdgEdgeType.Clamped, NdgEdgeType.Clamped, ...
    NdgEdgeType.Clamped, NdgEdgeType.Clamped];
            
if (type == NdgCellType.Tri)
    mesh = makeUniformTriUMeshUnion(N, [-1, 1], [-1, 1], ...
        M, M, bctype);
elseif(type == NdgCellType.Quad)
    mesh = makeUniformQuadUMeshUnion(N, [-1, 1], [-1, 1], ...
        M, M, bctype);
else
    msgID = 'AdvUniformUnion:inputCellTypeError';
    msgtext = ['The input cell type should be NdgCellType.Tri',...
        ' or NdgCellType.Tri.'];
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

