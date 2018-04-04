%> Khan AA, Lai W. Modeling shallow water flows using the discontinuous
%> Galerkin method. CRC Press; 2014.
%> 
classdef ParabolicBowl1d < SWEConventional1d
    
    properties(Constant)
        hmin = 1e-3
        gra = 9.81
        a = 3000
        h0 = 10
        B = 5
    end
    
    methods
        function obj = ParabolicBowl1d(N, M)
            obj = obj@SWEConventional1d();
            [ mesh ] = makeUniformMesh( N, M );
            obj.initPhysFromOptions( mesh );
            obj.fext = obj.evaluateExactFunc( obj.getOption('finalTime') );
        end
        
        function drawDepthError( obj )
            ncfile = [obj.getOption('outputNetcdfCaseName'), '.1-1.nc'];
            time = ncread( ncfile, 'time' );
            pos = makeNdgPostProcessFromNdgPhys( obj );
            err = zeros( pos.Nt, 1 );
            for t = 1:numel( time )
                fext = obj.evaluateExactFunc( time(t) );
                fphys = pos.accessOutputResultAtStepNum( t );
                temp = pos.evaluateNormErr2( fphys, fext );
                err( t ) = temp(1);
            end
            figure(1);
            plot( time, err, 'LineWidth', 2 );
            grid on; 
            hold on;
            xlabel('time (s)', 'Interpreter', 'latex','FontSize', 16);
            ylabel('$L^2(h)$', 'Interpreter', 'latex','FontSize', 16);
        end
        
        function drawFlowError( obj )
            ncfile = [obj.getOption('outputNetcdfCaseName'), '.1-1.nc'];
            time = ncread( ncfile, 'time' );
            pos = makeNdgPostProcessFromNdgPhys( obj );
            err = zeros( pos.Nt, 1 );
            for t = 1:numel( time )
                fext = obj.evaluateExactFunc( time(t) );
                fphys = pos.accessOutputResultAtStepNum( t );
                temp = pos.evaluateNormErr2( fphys, fext );
                err( t ) = temp(2);
            end
            figure(2);
            plot( time, err./1.1, 'LineWidth', 2 );
            grid on; 
            hold on;
            xlabel('time (s)', 'Interpreter', 'latex','FontSize', 16);
            ylabel('$L^2(hu)$', 'Interpreter', 'latex','FontSize', 16);
        end
        
        function stime = drawSection( obj, periodFrac )
            ncfile = [obj.getOption('outputNetcdfCaseName'), '.1-1.nc'];
            time = ncread( ncfile, 'time' );
            pos = makeNdgPostProcessFromNdgPhys( obj );

            stime = zeros( numel(periodFrac), 1 );
            for i = 1:numel(periodFrac)
                extTime = periodFrac(i) * obj.getOption('finalTime');
                [~, t] = min( abs(time - extTime) );
                stime(i) = time(t);
                fext = obj.evaluateExactFunc( time(t) );
                fphys = pos.accessOutputResultAtStepNum( t );
                
                figure( 2*(i-1) + 1 );
                drawElavationSection( obj.meshUnion(1).x, ...
                    fphys{1}, fext{1}, obj.fphys{1}(:,:,3) );
                figure( 2*i ); hold on;
                drawFluxSection( obj.meshUnion(1).x, ...
                    fphys{1}, fext{1} );
            end
        end
        
        function drawMassError( obj )
            pos = makeNdgPostProcessFromNdgPhys( obj );
            pos.checkMassVolume( 1 );
        end
        
        function drawSectionError( obj, timefrac )
            ncfile = [obj.getOption('outputNetcdfCaseName'), '.1-1.nc'];
            time = ncread( ncfile, 'time' );
            pos = makeNdgPostProcessFromNdgPhys( obj );
            stime = zeros( numel(timefrac), 1 );
            
            for i = 1:numel(timefrac)
                extTime = timefrac(i) * obj.getOption('finalTime');
                [~, t] = min( abs(time - extTime) );
                stime(i) = time(t);
                fext = obj.evaluateExactFunc( time(t) );
                fphys = pos.accessOutputResultAtStepNum( t );
                
                for m = 1:obj.Nmesh
                    mesh = obj.meshUnion(m);
                    temp = fphys{m}(:, :, 1) - fext{m}(:, :, 1);
                    figure(1)
                    plot( mesh.x(:), temp(:)); hold on;
                    
                    temp = fphys{m}(:, :, 2) - fext{m}(:, :, 2);
                    figure(2)
                    plot( mesh.x(:), temp(:)); hold on;
                    
                end
            end
        end
    end
    
    methods( Access = protected )
        
        function fphys = setInitialField( obj )
            fphys = obj.evaluateExactFunc( 0 );
        end
        
        function [ option ] = setOption( obj, option )
            T = 2*pi*obj.a/sqrt(2*obj.gra*obj.h0);
            ftime = 2*T;
            outputIntervalNum = 2000;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('cfl') = 1/obj.meshUnion(1).cell.N;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK33;
            option('limiterType') = NdgLimiterType.TVB;
            option('limiterParameter') = 1e-6;
            option('SWELimiterType') = SWELimiterType.OnDepth;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
        end
        
        function fext = evaluateExactFunc( obj, time )
            fext = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                
                fext{m} = zeros(mesh.cell.Np, mesh.K, obj.Nfield);
                % parameters
                g = obj.gra;            
                % Init Condition
                fext{m}(:,:,3) = obj.h0.*(mesh.x.^2./obj.a^2 - 1);
                w = sqrt(2*g*obj.h0)./obj.a;
                z = (- obj.B.^2 * (cos(2*w*time) + 1) ...
                    - (4*obj.B*w).* mesh.x * cos(w*time))./(4*g);
                fext{m}(:,:,1) = max( 0, z - fext{m}(:,:,3));
                fext{m}(:,:,2) = fext{m}(:,:,1).* obj.B.* sin(w*time);
            end
        end
    end
    
end

function drawFluxSection( x, fphys, fext )
ind = fphys( :, :, 1 ) < 1e-2;
u = fphys( :, :, 2 )./fphys(:, :, 1); u( ind ) = 0;
ue = fext( :, :, 2 )./fext(:, :, 1);  ue( ind ) = 0;
hold on;
p1 = plot( x(:), u(:), 'b.-', 'LineWidth', 2, 'MarkerSize', 4 );
p2 = plot( x(:), ue(:), 'r--', 'LineWidth', 2, 'MarkerSize', 4  );
% p1 = plot( x, fphys(:, 2), 'b', 'LineWidth', 1.5 );
% p2 = plot( x, fext(:, 2), 'r--', 'LineWidth', 1.5 );
grid on;
box on;
ylim([ -50, 50 ]);
xlabel('$x$ (m)', 'FontSize', 16, 'Interpreter', 'Latex');
ylabel('$u$ (m/s)', 'FontSize', 16, 'Interpreter', 'Latex');
% legend([p1, p2], {'Quadrature-free DG', 'Exact'}, ...
%     'FontSize', 16, 'Interpreter', 'Latex', 'box', 'off');
end

function drawElavationSection( x, fphys, fext, bot )
plot( x(:), bot(:), 'k' );
temp1 = fphys(:, :, 1) + bot;
temp2 = fext(:, :, 1) + bot;
hold on;
p1 = plot( x(:), temp1(:), 'b.-', 'LineWidth', 2, 'MarkerSize', 4  );
p2 = plot( x(:), temp2(:), 'r--', 'LineWidth', 2, 'MarkerSize', 4  );
grid on;
box on;
xlabel('$x$ (m)', 'FontSize', 16, 'Interpreter', 'Latex');
ylabel('$\eta$ (m)', 'FontSize', 16, 'Interpreter', 'Latex');
% legend([p1, p2], {'Quadrature-free DG', 'Exact'}, ...
%     'FontSize', 16, 'Interpreter', 'Latex', 'box', 'off');
end

function [ mesh ] = makeUniformMesh( N, M )
xlim = [-5000, 5000];
bcType = [NdgEdgeType.ZeroGrad, NdgEdgeType.ZeroGrad];
[ mesh ] = makeUniformMesh1d( N, xlim, M, bcType );
end
