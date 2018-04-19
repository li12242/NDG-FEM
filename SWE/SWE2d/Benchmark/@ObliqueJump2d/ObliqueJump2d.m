classdef ObliqueJump2d < SWEConventional2d
    properties (Constant)
        hmin = 1e-4
        h0 = 1.0
        u0 = 8.57
        gra = 9.81
    end
    
    methods ( Access = public )
        function obj = ObliqueJump2d( N )
            [ path, ~, ~ ] = fileparts( mfilename('fullpath') );
            meshfile = [ path, '/mesh/tri.msh' ];
            mesh = makeGmshFileUMeshUnion2d( N, meshfile );
            obj.initPhysFromOptions( mesh );
            obj.getExactFunc(  );
        end
    end% methods

    methods (Access = protected)
        function setOpenBoundaryCondition(obj)
            for m = 1:obj.Nmesh
                obj.fext{m}(:, :, 1) = obj.h0;
                obj.fext{m}(:, :, 2) = obj.h0 .* obj.u0;
            end
        end

        function [ option ] = setOption( obj, option )
            ftime = 20;
            outputIntervalNum = 200;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK22;
            option('limiterType') = NdgLimiterType.Vert;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('FrictionType') = SWEFrictionType.None;
        end

        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:, :, 1) = obj.h0;
                fphys{m}(:, :, 2) = obj.h0 .* obj.u0;
            end
        end
    end

    methods ( Access = private )
        %> get the exact field into fext variable
        function getExactFunc( obj )
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                obj.fext{m}(:, :, 1) = 1.0;
                k = -tan(30/180*pi);
                ind = ( mesh.y >= ( 30+k*( mesh.x-10 ) ));
                obj.fext{m}( ind ) = 1.5049;
            end
        end% func
    end
end