classdef AdvRotationUniformMesh2d < AdvAbstractVarFlow2d
    
    properties
        %> Number of basis function
        N
        %> Number of elements on each axis
        M
    end
    
    properties(Constant)
        %> domain central
        x0 = 0
        %> domain central
        y0 = 0
        %> distance of the initial Gauss mount from the central point
        rd = 0.5
        %> size of the initial Gauss mount
        r0 = 0.25
        %> angle velocity
        w = 5*pi/6;
    end
    
    methods
        function obj = AdvRotationUniformMesh2d(N, M, cellType)
            mesh = makeUniformMesh(N, M, cellType);
            obj = obj@AdvAbstractVarFlow2d( );
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
            option('finalTime') = 1.2;
            option('outputType') = enumOutputFile.NetCDF;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = 2.4/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('timeInterval') = sqrt(2)/obj.M/obj.w/(2*obj.N + 1);
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('limiterType') = enumLimiter.None;
        end
        
        %> the exact function
        function f_ext = getExtFunc(obj, x, y, time)
            [ Np, K ] = size( x );
            f_ext = zeros( Np, K, obj.Nfield );
            
            theta0 = -pi;
            theta = theta0 + obj.w * time;
            xt = obj.x0 + obj.rd*cos(theta);
            yt = obj.y0 + obj.rd*sin(theta);
            r2 = sqrt( (x - xt).^2+(y - yt).^2 )./obj.r0;
            ind = ( r2 <= 1.0);
            temp = zeros( Np, K );
            temp(ind) = ( 1+cos( r2(ind)*pi ) )./2;
            
            f_ext(:,:,1) = temp;
            f_ext(:,:,2) = obj.w.* (- y);
            f_ext(:,:,3) = obj.w.*( x );
        end% func
    end% methods
    
end% class

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

