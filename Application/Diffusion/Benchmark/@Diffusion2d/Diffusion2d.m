classdef Diffusion2d < DiffusionAbstract2d
    %DIFFUSION2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties ( Constant )
        x0 = 0
        y0 = 0
        miu = 1e-1;
    end
    
    methods
        function obj = Diffusion2d( N, M, cellType )
            mesh = makeUniformMesh(N, M, cellType);
            obj = obj@DiffusionAbstract2d( );
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
        
        function option = setOption( obj, option )
            outputIntervalNum = 50;
            option('startTime') = 0.0;
            option('finalTime') = 0.5;
            option('outputType') = enumOutputFile.NetCDF;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = 2.4/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('timeInterval') = (2 / obj.M / (obj.N+1)).^2 ./ obj.miu / 100;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('limiterType') = enumLimiter.None;
        end
        
        function f = getExtFunc( obj, x, y, time )
            % f = x + y;
            f = exp( - ( x.^2 + y.^2 )./ obj.miu / (4 * time + 1) ) / ( 4 * time + 1 );
        end
    end
    
end

function mesh = makeUniformMesh(N, M, type)
bctype = [...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped];

if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, [-1, 1], [-1, 1], M, M, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N, [-1, 1], [-1, 1], M, M, bctype);
else
    msgID = [mfilename, ':inputCellTypeError'];
    msgtext = ['The input cell type should be NdgCellType.Tri',...
        ' or NdgCellType.Tri.'];
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func


