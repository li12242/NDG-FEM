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
                mesh = obj.meshUnion(m);
                fphys{m} = getExtFunc(obj, mesh.x, mesh.y, 0);
            end
        end% func
        
        function [ option ] = setOption( obj, option )
            outputIntervalNum = 50;
            finalTime = 2.0;
            option('startTime') = 0.0;
            option('finalTime') = finalTime;
            option('timeInterval') ...
                = 2/obj.M/sqrt(obj.u0 ^ 2 + obj.v0 ^2)/(2*obj.N + 1);
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = finalTime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
        end
               
        %> the exact function
        function [ f_ext ] = getExtFunc(obj, x, y, time)
            xc = obj.x0 + obj.u0.*time;
            yc = obj.y0 + obj.v0.*time;
            
            sigma = 125*1e3/(33*33);
            t = -( (x-xc).^2+(y-yc).^2 )*sigma;
            f_ext = exp(t);
        end% func
    end
    
end

function mesh = makeUniformMesh(N, M, type)
bctype = [enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped];

if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, [-1, 1], [-1, 1], ...
        M, M, bctype);
elseif(type == enumStdCell.Quad)
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

