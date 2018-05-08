%> The advection problem with the contant flow filed.
%> 
classdef ConstAdvUniformMesh2d < AdvAbstractConstFlow2d
    
    properties
        %> Number of basis function
        N
        %> Number of elements on each axis
        M
    end
    
    properties( Constant )
        x0 = -0.5;
        y0 = -0.5;
        u0 = 0.5;
        v0 = 0.5;
    end
    
    methods
        function obj = ConstAdvUniformMesh2d(N, M, cellType)
            mesh = makeUniformMesh(N, M, cellType);
            obj = obj@AdvAbstractConstFlow2d();
            obj.N = N;
            obj.M = M;
            obj.initPhysFromOptions( mesh );            
        end
    end
    
    methods( Access = protected )
        
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                fphys{m} = getExtFunc(obj, obj.meshUnion(m), 0);
            end
        end% func
        
        function [ option ] = setOption( obj, option )
            outputIntervalNum = 50;
            finalTime = 2.0;
            option('startTime') = 0.0;
            option('finalTime') = finalTime;
            option('temporalDiscreteType') = NdgTemporalIntervalType.Constant;
%             option('CFL') = 0.2;
            option('timeInterval') ...
                = 2/obj.M/sqrt(obj.u0 ^ 2 + obj.v0 ^2)/(2*obj.N + 1);
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = finalTime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('limiterType') = NdgLimiterType.None;
        end
               
        %> the exact function
        function [ f_ext ] = getExtFunc(obj, mesh, time)
            xc = obj.x0 + obj.u0.*time;
            yc = obj.y0 + obj.v0.*time;
            
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
    mesh = makeUniformTriMesh(N, [-1, 1], [-1, 1], ...
        M, M, bctype);
elseif(type == NdgCellType.Quad)
    mesh = makeUniformQuadMesh(N, [-1, 1], [-1, 1], ...
        M, M, bctype);
else
    msgID = [mfilename, ':inputCellTypeError'];
    msgtext = ['The input cell type should be NdgCellType.Tri',...
        ' or NdgCellType.Tri.'];
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

