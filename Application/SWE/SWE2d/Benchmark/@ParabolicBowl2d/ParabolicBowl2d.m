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
        mass
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
        
        function massError( obj )
            pos = makeNdgPostProcessFromNdgPhys( obj );
            pos.checkMassVolume( 1 );
        end
        
        function depthError( obj )
            ncfile = [obj.getOption('outputNetcdfCaseName'), '.1-1.nc'];
            time = ncread( ncfile, 'time' );
            pos = makeNdgPostProcessFromNdgPhys( obj );
            err = zeros( pos.Nt, 1 );
            for t = 1:numel( time )
                fext = obj.getExactFunction( time(t) );
                fphys = pos.accessOutputResultAtStepNum( t );
                temp = pos.evaluateNormErr2( fphys, fext );
                err( t ) = temp(3);
            end
            plot( time, err );
        end
        
        function stime = drawSection( obj, periodFrac )
            ncfile = [obj.getOption('outputNetcdfCaseName'), '.1-1.nc'];
            time = ncread( ncfile, 'time' );
            pos = makeNdgPostProcessFromNdgPhys( obj );
            Ng = 50;
            xg = linspace( -4e3, 4e3, Ng )';
            yg = zeros( Ng, 1 );
            zg = yg;
            bot = obj.a.*( xg.^2 + yg.^2 );
            stime = zeros( numel(periodFrac), 1 );
            for i = 1:numel(periodFrac)
                extTime = periodFrac(i) * obj.T;
                [~, t] = min( abs(time - extTime) );
                stime(i) = time(t)./obj.T;
                fg = pos.interpolateOutputStepResultToGaugePoint( xg, yg, zg, t );
                fext = obj.getExactFunction( time(t) );
                feg = pos.interpolatePhysFieldToGaugePoint( fext, xg, yg, zg );
                drawElavationSection( xg, fg, feg, bot );
                drawFluxSection( xg, fg, feg );
            end
        end
        
    end
    
    methods(Access=protected)
        function fphys = setInitialField( obj )
            fphys = obj.getExactFunction(0);
        end
        
        function fext = setExtField( obj )
            fext = getExactFunction(obj, obj.getOption('finalTime') );
        end
        
        function [ option ] = setOption( obj, option )
            obj.T = 2*pi./sqrt(8*obj.gra*obj.a);
            ftime = obj.T;
            
            outputIntervalNum = 200;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
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

function drawFluxSection( x, fphys, fext )
figure('Color', 'w'); hold on;
u = zeros( size(fphys, 1), 1 );
ue = zeros( size(fext, 1), 1 );
ind = (fphys(:, 1) > 1e-4); u(ind) = fphys(ind, 2)./fphys(ind, 1);
ind = (fphys(:, 1) > 1e-4); ue(ind) = fext(ind, 2)./fext(ind, 1);
p1 = plot( x, u, 'b', 'LineWidth', 1.5  );
p2 = plot( x, ue, 'r--', 'LineWidth', 1.5  );
% p1 = plot( x, fphys(:, 2), 'b', 'LineWidth', 1.5 );
% p2 = plot( x, fext(:, 2), 'r--', 'LineWidth', 1.5 );
grid on;
box on;
ylim([ -4, 4 ]);
xlabel('$x$(m)', 'FontSize', 16, 'Interpreter', 'Latex');
ylabel('$u$(m/s)', 'FontSize', 16, 'Interpreter', 'Latex');
legend([p1, p2], {'Quadrature-free DG', 'Exact'}, ...
    'FontSize', 16, 'Interpreter', 'Latex', 'box', 'off');
end

function drawElavationSection( x, fphys, fext, bot )
figure('Color', 'w'); hold on;
plot( x, bot, 'k' );
p1 = plot( x, fphys(:, 1) + bot, 'b', 'LineWidth', 1.5  );
p2 = plot( x, fext(:, 1) + bot, 'r--', 'LineWidth', 1.5  );
grid on;
box on;
xlabel('$x$(m)', 'FontSize', 16, 'Interpreter', 'Latex');
ylabel('$\eta$(m)', 'FontSize', 16, 'Interpreter', 'Latex');
legend([p1, p2], {'Quadrature-free DG', 'Exact'}, ...
    'FontSize', 16, 'Interpreter', 'Latex', 'box', 'off');
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

