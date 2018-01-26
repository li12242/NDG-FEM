classdef TriHump1d < SWEWD1d
    
    properties(Constant)
        h0 = 0.75
        n = (1.25e-2)^2 % Manning
        hmin = 2e-3
        dam_pos = 15.5;
        gra = 9.81;
    end
    
    methods
        function obj = TriHump1d( N, M )
            obj = obj@SWEWD1d();
            [ mesh ] = makeUniformMesh( N, M );
            obj.initPhysFromOptions( mesh );
        end
        
        function drawGaugeDepth( obj )
            name = {'p1.txt', 'p2.txt', 'p3.txt', ...
                'p4.txt', 'p5.txt', 'p6.txt', 'p7.txt'};
            [path, ~, ~] = fileparts(mfilename('fullpath'));
            humpos = makeNdgPostProcessFromNdgPhys(obj);
            xg = 15.5 + [2, 4, 8, 10, 11, 13, 20];
            [ gaugeValue ] = humpos.interpolateOutputResultToGaugePoint( xg, xg, xg );
            time = ncread( humpos.outputFile{1}, 'time');
            for i = 1:numel( name )
                filename = fullfile(path, name{i});
                [ t, h ] = read_measured_data(filename);
                figure(i);
                plot(t, h, 'ro'); hold on;
                temp = gaugeValue(i, 1, :);
                plot(time, temp(:), 'g-')
                xlim([0, 90]); ylim([0, .75]);
                grid on; box on;
                xlabel('time (s)', 'Interpreter', 'Latex', 'FontSize', 16);
                ylabel('h (m)', 'Interpreter', 'Latex', 'FontSize', 16);
                legend({'Measured', 'Numerical Result'}, 'box', 'off', 'FontSize', 16, ...
                    'Interpreter', 'Latex')
            end
            
        end
    end
    
    methods( Access = protected )
        
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                % initial condition
                fphys{m} = zeros(mesh.cell.Np, mesh.K, obj.Nfield);
                fphys{m}(:, mesh.xc < obj.dam_pos, 1) = obj.h0;
                % bottom topography
                bot = zeros( mesh.cell.Np, mesh.K );
                x0 = 28.5; len = 6/2; b0 = 0.4;
                ind = ( mesh.x > (x0-len) ) & ( mesh.x < (x0+len) );
                bot(ind) = b0 - abs( mesh.x(ind) - x0 )*b0/len;
                fphys{m}(:, :, 3) = bot;
            end
        end
        
        function [ option ] = setOption( obj, option )            
            ftime = 93;
            outputIntervalNum = 80;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('cfl') = 1/obj.meshUnion(1).cell.N;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.TVB;
            option('limiterParameter') = 1e-10;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
        end
    end
    
end

function [ t, h ] = read_measured_data( filename )
fp = fopen(filename);
for i = 1:3
    fgetl(fp);
end
Np = fscanf(fp, '%d', 1);
data = fscanf(fp, '%f,%f\n', [2, Np]);
fclose(fp);

t = data(1, :);
h = data(2, :);
end

function [ mesh ] = makeUniformMesh( N, M )
xlim = [0, 38];
bcType = [NdgEdgeType.SlipWall, NdgEdgeType.ZeroGrad];
[ mesh ] = makeUniformMesh1d( N, xlim, M, bcType );
end
