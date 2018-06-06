classdef PartialDamBreak < SWEConventional2d
    properties (Constant)
        hmin = 1e-4
        h0 = 5.0
        h1 = 10
        gra = 9.81
    end

    methods ( Access = public )
        function obj = PartialDamBreak( N )
            [ path, ~, ~ ] = fileparts( mfilename('fullpath') );
            meshfile = [ path, '/mesh/quad.msh' ];
            mesh = makeGmshFileUMeshUnion2d( N, meshfile );
            obj.initPhysFromOptions( mesh );
        end
    end% methods

    methods (Access = protected)
        function [ option ] = setOption( obj, option )
            ftime = 7;
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
        end% func

        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                ind = (mesh.EToR == NdgRegionType.Wet);
                fphys{m}(:, :, 1) = 0;
                fphys{m}(:, ind, 1) = obj.h1;
                
            end
        end% func
    end% methods
end% class