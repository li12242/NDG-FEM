classdef AdvRotationUniformMesh3d < AdvAbstractVarFlow3d
    %ADVROTATIONUNIFORMMESH3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %> max order of horizontal and vertical basis functions
        N, Nz
        %> num of elements on horizontal axis
        M, Mz
    end
    
    properties (Constant)
        %> domain central
        x0 = 0, 
        y0 = 0
        %> distance of the initial Gauss mount from the central point
        rd = 0.5
        %> size of the initial Gauss mount
        r0 = 0.25
        %> angle velocity
        w = 5 * pi / 6;
    end
    
    methods
        function obj = AdvRotationUniformMesh3d(N, Nz, M, Mz, cellType)
            mesh = makeUniformMesh(N, Nz, M, Mz, cellType);
            obj = obj@AdvAbstractVarFlow3d( );
            obj.N = N;
            obj.M = M;
            obj.Nz = Nz;
            obj.Mz = Mz;
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
        
        function option = setOption( obj, option )
            outputIntervalNum = 100;
            option('startTime') = 0.0;
            option('finalTime') = 2;
            option('outputType') = enumOutputFile.NetCDF;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = 2.4/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            dtz = 2/obj.Mz/0.5/(2*obj.Nz + 1);
            dth = sqrt(2)/obj.M/obj.w/(2*obj.N + 1);
            option('timeInterval') = min( dtz, dth );
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('limiterType') = enumLimiter.None;
        end
        
        %> the exact function
        function f_ext = getExtFunc( obj, mesh, time )
            f_ext = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
            
            theta0 = - pi;
            theta = theta0 + obj.w * time;
            xt = obj.x0 + obj.rd * cos(theta);
            yt = obj.y0 + obj.rd * sin(theta);
            zt = 0.5;
            r2 = sqrt((mesh.x - xt).^2 + (mesh.y - yt).^2 + (mesh.z - zt).^2 ) ./ obj.r0;
            ind = ( r2 <= 1.0);
            temp = zeros( mesh.cell.Np, mesh.K );
            temp(ind) = ( 1 + cos( r2(ind) * pi ) ) ./ 2;
            
            f_ext(:,:,1) = temp;
            f_ext(:,:,2) = obj.w .* ( - mesh.y );
            f_ext(:,:,3) = obj.w .* ( + mesh.x );
            f_ext(:,:,4) = - 0.5;
        end% func
    end% methods
end

function mesh = makeUniformMesh(N, Nz, M, Mz, type)
    bctype = [...
        enumBoundaryCondition.Clamped, ...
        enumBoundaryCondition.Clamped, ...
        enumBoundaryCondition.Clamped, ...
        enumBoundaryCondition.Clamped];

    if (type == enumStdCell.PrismTri)
        mesh2d = makeUniformTriMesh(N, [-1, 1], [-1, 1], M, M, bctype);
        std = StdPrismTri( N, Nz );
    elseif(type == enumStdCell.PrismQuad)
        mesh2d = makeUniformQuadMesh(N, [-1, 1], [-1, 1], M, M, bctype);
        std = StdPrismQuad( N, Nz );
    end

    zs = ones(mesh2d.Nv, 1);
    zb = ones(mesh2d.Nv, 1) - 2.0;
    mesh = NdgExtendMesh3d( std, mesh2d, zs, zb, Mz );
    mesh.InnerSideEdge = NdgSideEdge3d( mesh, 1 );
    mesh.BottomEdge = NdgBottomEdge3d( mesh, 1 );
end% func

