classdef ConstAdvGmshMesh2d < AdvAbstractConstFlow2d
    %CONSTADVGMSHMESH2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties( Constant )
        x0 = -0.5;
        y0 = -0.5;
        u0 = 0.5;
        v0 = 0.5;
    end
    
    properties
        N
        gmshFile
    end
    
    methods
        function obj = ConstAdvGmshMesh2d( N )
            obj = obj@AdvAbstractConstFlow2d();
            obj.N = N;
            obj.gmshFile = [fileparts( mfilename('fullpath') ), '/mesh/TriMesh.msh'];
            mesh = makeGmshFileUMeshUnion2d( N, obj.gmshFile );
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
            option('startTime') = 0.0;
            option('finalTime') = 2.0;
            option('temporalDiscreteType') = NdgTemporalIntervalType.Constant;

            dx = 2;
            for m = 1:obj.Nmesh
                dx = min( dx, min( obj.meshUnion(m).charLength ) );
            end
            option('timeInterval') = dx/sqrt(obj.u0.^2 + obj.v0.^2)/(2*obj.N + 1);
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = 2/outputIntervalNum;
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
    

