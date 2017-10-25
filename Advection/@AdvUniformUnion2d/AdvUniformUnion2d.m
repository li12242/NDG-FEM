classdef AdvUniformUnion2d < Adv2d
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
        %> final time
        finalTime = 2
        %> start time
        initialTime = 0
        
    end
    
    methods
        
        function obj = AdvUniformUnion2d(N, M, type)
            % check input and create uniform mesh
            mesh = makeUniformMesh(N, M, type);
            solver = NdgSolver();
            obj = obj@Adv2d(solver, mesh);
            obj.initSolver(solver, mesh);
            obj.initPhysField;
        end% func
        
        function initSolver(obj, solver, mesh)
            udet = sqrt( obj.u0.^2 + obj.v0^2 );
            dx = max( mesh.x(:, 1) ) - min( mesh.x(:, 1) );
            dt = dx/udet/(2*mesh.cell.N + 1);
            
            solver.setTimeIntervalType( NdgIntervalType.Constant,  dt);
            solver.setLimiterType( NdgLimiterType.None );
            solver.setTemporalDiscreteType( NdgTemporalDiscreteType.RK45 );
            solver.setBoundaryType( NdgBCType.None );
            solver.setOutputFile( 'advUniformMesh', ...
                NdgIntervalType.DeltaTime, 0.002 );
        end
        
        function initPhysField(obj)
            mesh = obj.mesh;
            Nmesh = obj.Nmesh;
            K = mesh.K;
            Np = mesh.cell.Np;
            obj.fphy = zeros(Np, K, obj.Nfield, Nmesh);
            obj.fphy(:, :, 1, 1) = obj.ext_func(0);
            obj.fphy(:, :, 2, 1) = obj.u0;
            obj.fphy(:, :, 3, 1) = obj.v0;
        end% func
    end
    
    methods(Access = protected)
        %> the exact function
        function f_ext = ext_func(obj, time)
            xc = obj.x0 + obj.u0.*time;
            yc = obj.y0 + obj.v0.*time;
            
            mesh = obj.mesh;
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

