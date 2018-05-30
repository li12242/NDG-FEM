classdef Sanya < SWEPreBlanaced2d

    
    properties( SetAccess = protected )
    end
    
    properties(Constant)
         %> wet/dry depth threshold
        hmin = 0.5        
        %> gravity acceleration
        gra = 9.8        
        %> interval of tide elevation
        tideinterval = 600
        
    end
    
    properties
        N
        gmshFile
        OBVid
        Tide
        
        FextLoader
    end
    
    methods 
        function obj = Sanya( N                                                                                                                                                                                                                                                                                      )
            obj = obj@SWEPreBlanaced2d();
            obj.N = N;
            obj.gmshFile = [pwd, '/SWE2d/'...
                '@Sanya/mesh/sanya0111'];
            mesh = makeGmshFileUMeshUnion2d( N, obj.gmshFile );
            obj.initPhysFromOptions( mesh );
            obj.FindOBVertex;
            obj.ReadTideElevation;
        end 
        
        
    end%methods
    
    methods(Access=protected)
        function fphys = setInitialField( obj )
            fphys = getInitialFunction(obj);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 259200;
            outputIntervalNum = 432;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('temporalDiscreteType') = NdgTemporalIntervalType.DeltaTime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = 'Sanya2k_0328';
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.Vert;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('CoriolisType')=CoriolisType.Latitude;
            option('LatitudeFilePath')='SWE2d/@Sanya/tide/vertex_lat.txt';
            option('WindType')=WindType.None;
            option('FrictionType')=FrictionType.Quad;
            option('FrictionCoefficient_n')=0.017;
%             option('ExternalFieldInterpolateType') = FextInterpolateType.Linear;
        end
            
%         matEvaluateRK45( obj );
        matUpdateExternalField( obj, time, fphys )
        
    end%methods
    
end%classdef