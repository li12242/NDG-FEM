classdef OpenChannel2d < SWEConventional2d
    
    properties(Constant)
        hmin = 1e-4; % water depth threshold
        gra = 9.8; % gravity acceleration
        %> channel depth
        H = 320;
        %> amplitude
        eta = 0.2;
        %> wave period
        T = 500;
        
        %> domain length
        ChLength = 200e3;
        %> domain width
        ChWidth = 500;
        %> mesh rotation angle
        theta = 0;
        %> mesh center
        xc = 0;
        %> mesh center
        yc = 0;
    end
    
    methods
        function obj = OpenChannel2d( N, M, type )
            obj = obj@SWEConventional2d();
            mesh = obj.makeUniformMesh( N, M, type );
            obj.initPhysFromOptions( mesh );
        end
        
        %> draw the numerical results at gauge points
        drawGaugePoints( obj, xg );
        
        function showParameter( obj )
            waveLen = obj.T * sqrt(obj.gra * obj.H);
            fprintf('\nWave parameter:\n');
            fprintf('\tWave amplitude (m): %f\n', obj.eta);
            fprintf('\tWave period (s): %f\n', obj.T);
            fprintf('\tWave speed (m/s): %f\n', sqrt(obj.gra * obj.H ));
            fprintf('\tWave length (m): %f\n', waveLen);
            
            fprintf('\nChannel parameter\n');
            fprintf('\tChannel depth (m): %f\n', obj.H);
            fprintf('\tChannle length (m): %f\n', obj.ChLength);
            fprintf('\tNumber of wave: %f\n', obj.ChLength / waveLen);
            fprintf('\tTime to propage through the channel (h): %f\n', ...
                obj.ChLength/sqrt(obj.gra * obj.H ) / 3600 );
            
            Ne = obj.meshUnion.K/2;
            fprintf('\nModel parameter\n');
            fprintf('\tElement number: %d\n', Ne );
            fprintf('\tNumber of elements in each wave: %f\n',  Ne/(obj.ChLength / waveLen));
            fprintf('\tSimulate time (h) = %f\n', obj.getOption('finalTime')/3600 );
            fprintf('\tNumber of waves pass through the channel: %f\n', ...
                obj.getOption('finalTime')/( obj.ChLength/sqrt(obj.gra * obj.H ) ) );
            fprintf('\tOutput interval: %d\n', obj.getOption('outputTimeInterval'));
        end
    end
    
    methods( Access = protected, Static  )
        %> set open boundary condition
        obtype = setOpenBoundaryCondition( )
    end
    
    methods( Access = protected )
        function mesh = makeUniformMesh( obj, N, M, type )
            obtype = obj.setOpenBoundaryCondition();
            bctype = [NdgEdgeType.SlipWall, NdgEdgeType.SlipWall, ...
                obtype(1), obtype(2)];
            
            if (type == NdgCellType.Tri)
                mesh = makeUniformTriMesh(N, [0, obj.ChLength], [0, obj.ChWidth], ...
                    M, 2, bctype);
            elseif(type == NdgCellType.Quad)
                mesh = makeUniformQuadMesh(N, [0, obj.ChLength], [0, obj.ChWidth], ...
                    M, 2, bctype);
            else
                msgID = [mfilename, ':inputCellTypeError'];
                msgtext = ['The input cell type should be NdgCellType.Tri',...
                    ' or NdgCellType.Tri.'];
                ME = MException(msgID, msgtext);
                throw(ME);
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 7200 * 2;
            outputIntervalNum = 2000;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('temporalDiscreteType') = NdgTemporalIntervalType.DeltaTime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = [mfilename, '.', num2str(obj.meshUnion.cell.N)];
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('limiterType') = NdgLimiterType.None;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('integralType') = NdgDiscreteIntegralType.GaussQuadrature;
            option('CoriolisType') = SWECoriolisType.None;
            option('WindType') = SWEWindType.None;
            option('FrictionType') = SWEFrictionType.None;
        end
        
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                [ fphys{m}(:,:,1), fphys{m}(:,:,2) ] = obj.setExactSolution( mesh.x, 0.0 );
                fphys{m}(:, :, 4) = - obj.H;
            end
        end
        
        function matUpdateExternalField( obj, time, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                [ obj.fext{m}(:,:,1), obj.fext{m}(:,:,2) ] ...
                    = obj.setExactSolution( mesh.x, time );
                fphys{m}(1, end, 1) = obj.H;
                fphys{m}(1, end, 2) = 0;
            end
        end
        
        function [h, hu] = setExactSolution( obj, x, time )
            w = 2 * pi / obj.T;
            c = sqrt( obj.gra * obj.H );
            k = w/c;
            
            temp = cos( k .* x - w * time );
            h = obj.eta * temp + obj.H;
            u = obj.eta * sqrt( obj.gra/obj.H ) * temp;
            % fphys(:, :, 1) = h;
            hu = h .* u;
        end% func
        
        
    end% methods
    
end

