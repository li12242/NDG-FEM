classdef FrictionParabolicBowl2d < SWEConventional2d
    properties (Constant)
        hmin = 1e-4
        gra = 9.81
    end

    properties (Constant)
        h0 = 10
        a = 3e3
        B = 5
        k = 0.002
        T = 2*672
    end

    methods ( Access = public )
        function obj = FrictionParabolicBowl2d(N, M, cellType)
            [ mesh ] = makeUniformMesh(N, M, cellType );
            obj = obj@SWEConventional2d();
            obj.initPhysFromOptions( mesh );
            
            obj.fext = obj.getExactFunction( obj.getOption('finalTime') );
        end
    end

    methods (Access = protected)
        function fphys = setInitialField( obj )
            fphys = obj.getExactFunction(0);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = obj.T;
            outputIntervalNum = 500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('temporalDiscreteType') = NdgTemporalIntervalType.DeltaTime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK22;
            option('limiterType') = NdgLimiterType.Vert;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
        end
    end

    methods (Access = private)
        function [ fext ] = getExactFunction( obj, time )
            fext = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fext{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                
                p = sqrt(8*obj.gra*obj.h0./obj.a^2);
                s = sqrt(p^2 - obj.k^2);

                eta = obj.h0 - obj.B^2./2./obj.gra*exp(-obj.k*time) ...
                    - obj.B./obj.gra*exp(-obj.k*time/2)*(obj.k/2*sin(s*time) + s*cos(s*time)).*mesh.x ...
                    - obj.B./obj.gra*exp(-obj.k*time/2)*(obj.k/2*cos(s*time) + s*sin(s*time)).*mesh.y;

                fext{m}(:, :, 4) = (mesh.x.^2 + mesh.y.^2).*(obj.h0/obj.a^2);
                
                h = eta - fext{m}(:, :, 4); h(h<0) = 0; 
                u = obj.B*exp(-obj.k*time/2)*sin(s*time);
                v = obj.B*exp(-obj.k*time/2)*sin(s*time);                

                fext{m}(:, :, 1) = h;
                fext{m}(:, :, 2) = h .* u;
                fext{m}(:, :, 3) = h .* v;
                
            end
        end
    end
end% class

function [ mesh ] = makeUniformMesh(N, M, type)
bctype = [...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad];

if (type == NdgCellType.Tri)
    mesh = makeUniformTriMesh(N, [-6000, 6000], [-6000, 6000], ...
        M, M, bctype );
elseif(type == NdgCellType.Quad)
    mesh = makeUniformQuadMesh(N, [-6000, 6000], [-6000, 6000], ...
        M, M, bctype );
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func