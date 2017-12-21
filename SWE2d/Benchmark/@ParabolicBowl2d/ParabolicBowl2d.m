%> Parabolic bowl bottom with oscillating flow.
%
%> This test case is to examine the capacity of the methods to treat
%> flooding and drying. Details of this example please refer to Ern et. al
%> (2008).
%>
%> [1]. Ern A, Piperno S, Djadel K. A well-balanced Runge-Kutta discontinuous
%> Galerkin method for the shallow-water equations with flooding and drying.
%> International Journal for Numerical Methods in Fluids 2008;58:1–25.
%> doi:10.1002/fld.1674.
%>
classdef ParabolicBowl2d < SWEConventional2d
    
    properties( Constant )
        hmin = 1e-4
        gra = 9.81
        a = 1.6e-7
        X = 1
        Y = -0.41884
    end
    
    properties(SetAccess=private)
        T
    end
    
    methods
        function obj = ParabolicBowl2d( varargin )
            if nargin == 1
                gmshFile = [pwd, ...
                    '/SWE/Benchmark/@ParabolicBowl2d/mesh/TriMesh.msh'];
                N = varargin{1};
                mesh = makeGmshFileUMeshUnion2d( N, gmshFile );
                
            elseif nargin == 3
                N = varargin{1};
                M = varargin{2};
                cellType = varargin{3};
                [ mesh ] = makeUniformMesh(N, M, cellType );
            end
            obj = obj@SWEConventional2d();
            obj.initPhysFromOptions( mesh );
            obj.fext = obj.setExtField(  );
        end
    end
    methods(Access=protected)
        function fphys = setInitialField( obj )
            fphys = getExactFunction(obj, 0);
        end
        
        function fext = setExtField( obj )
            fext = getExactFunction(obj, obj.getOption('finalTime') );
        end
        
        function [ option ] = setOption( obj, option )
            obj.T = 2*pi./sqrt(8*obj.gra*obj.a);
            ftime = obj.T/4;
            
            outputIntervalNum = 80;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('temporalDiscreteType') = NdgTemporalIntervalType.DeltaTime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.Vert;
            option('equationType') = NdgDiscreteEquationType.Weak;
            option('integralType') = NdgDiscreteIntegralType.GaussQuadrature;
        end
        
        function fphys = getExactFunction( obj, time )
            w  = sqrt(8*obj.gra.*obj.a);
            temp = obj.X+obj.Y*cos(w*time);
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                
                mesh = obj.meshUnion(m);
                fphys{m} = zeros(mesh.cell.Np, mesh.K, obj.Nfield);
                r2 = mesh.x.^2 + mesh.y.^2;
                h = 1./temp + obj.a*(obj.Y^2 - obj.X^2)*r2./temp.^2;
                h(h<0) = 0;
                fphys{m}(:,:,1) = h;
                
                u = - obj.Y*w*sin(w*time)./temp.*mesh.x./2;
                v = - obj.Y*w*sin(w*time)./temp.*mesh.y./2;
                
                fphys{m}(:, :, 2) = u.*h;
                fphys{m}(:, :, 3) = v.*h;
                fphys{m}(:, :, 4) = obj.a.*r2;
            end
        end
    end
    
end

function [ mesh ] = makeUniformMesh(N, M, type)
bctype = [...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad, ...
    NdgEdgeType.ZeroGrad];

if (type == NdgCellType.Tri)
    mesh = makeUniformTriMesh(N, [-4000, 4000], [-4000, 4000], ...
        M, M, bctype );
elseif(type == NdgCellType.Quad)
    mesh = makeUniformQuadMesh(N, [-4000, 4000], [-4000, 4000], ...
        M, M, bctype );
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

