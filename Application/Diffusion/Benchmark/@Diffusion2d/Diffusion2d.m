classdef Diffusion2d < DiffusionAbstract2d
    %DIFFUSION2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties ( Constant )
        len = 4
        x0 = 0
        y0 = 0
        miu = 1e-1;
    end
    
    methods
        function obj = Diffusion2d( N, M, cellType )
            obj = obj@DiffusionAbstract2d( );
            obj.N = N;
            obj.M = M;
            mesh = makeUniformMesh( obj, N, M, cellType);
            obj.initPhysFromOptions( mesh );
            
            Introduction(obj);
        end
    end
    
    methods( Access = protected )
        function Introduction( obj )
            fprintf('\n================  Diffusion2d  ==================');
            fprintf('\nDiffusion problem with constant viscosity.');
            fprintf('\nParameter: ');
            fprintf( '\n   casename = %s', obj.getOption('outputCaseName'));
            fprintf( '\n   miu = %8.4f', obj.miu);
            fprintf( '\n   dt  = %8.4f', obj.getOption('timeInterval'));
            fprintf( '\n   final time = %8.4f', obj.getOption('finalTime'));
            fprintf( '\n=================================================');
        end
        
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
            option('timeInterval') = (obj.len / obj.M / (obj.N+1)).^2 ./ obj.miu / 4;
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

function mesh = makeUniformMesh(obj, N, M, type)
bctype = [...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped];

domain = obj.len/2 * [-1, 1];
if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, domain, domain, M, M, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N, domain, domain, M, M, bctype);
else
    msgID = [mfilename, ':inputCellTypeError'];
    msgtext = ['The input cell type should be NdgCellType.Tri',...
        ' or NdgCellType.Tri.'];
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func


