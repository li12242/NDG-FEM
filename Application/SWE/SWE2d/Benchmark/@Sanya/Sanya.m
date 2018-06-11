classdef Sanya < SWEPreBlanaced2d
    
    properties(Constant)
        %> wet/dry depth threshold
        hmin = 0.5
        %> gravity acceleration
        gra = 9.8
        %> interval of tide elevation (s)
        tideinterval = 600
        %> tidal elevation of open boundary vertices
        tidalFile = ...
            'Application/SWE/SWE2d/Benchmark/@Sanya/tide/TideElevation.txt'
    end
    
    properties
        N
        %> open boundary edge index
        OBEdgeIndex
        %> open boundary tidal elevation
        Tide
    end
    
    methods
        function obj = Sanya( N                                                                                                                                                                                                                                                                                      )
            obj = obj@SWEPreBlanaced2d();
            obj.N = N;
            gmshFile = [ fileparts( mfilename('fullpath') ), '/mesh/sanya0111.msh' ];
            mesh = makeGmshFileUMeshUnion2d( N, gmshFile );
            obj.initPhysFromOptions( mesh );
        end
        
        OutputOpenBoundaryVertCoor( obj );
        ReadTideElevation( obj );
        
    end%methods
    
    methods ( Access = protected )
        function fphys = setInitialField( obj )
            fphys = getInitialFunction( obj );
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 259200;
            outputIntervalNum = 500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('temporalDiscreteType') = enumTemproalInterval.DeltaTime;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = 'Sanya2k_0328';
            option('limiterType') = enumLimiter.Vert;
            option('SWELimiterType') = enumSWELimiter.OnElevation;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('CoriolisType')= enumSWECoriolis.None;
            option('LatitudeFilePath')= ...
                [ fileparts( mfilename('fullpath') ),'/tide/vertex_lat.txt'];
            option('WindType')= enumSWEWind.None;
            option('FrictionType')= enumSWEFriction.None;
            option('FrictionCoefficient_n') = 0.017;
        end
        matUpdateExternalField( obj, time, fphys )
        
    end%methods
    
end%classdef